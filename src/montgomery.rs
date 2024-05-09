use std::fmt;
use std::fmt::Debug;

use crypto_bigint::{subtle::Choice, Encoding, Random, Zero, U512};
use rand::thread_rng;
use rug::Integer;
use rug::integer::Order;

use crate::constants::{ModP, HASSE_INTERVAL, PRIMES, PRIMES16, P_GMP};

#[derive(Clone, Debug, PartialEq)]
pub struct Point {
    pub x: ModP,
    pub z: ModP,
}

impl Point {
    pub const fn from_x(x: ModP) -> Point {
        Point { x, z: ModP::ONE }
    }

    pub const fn zero() -> Point {
        Point {
            x: ModP::ONE,
            z: ModP::ZERO,
        }
    }

    pub fn random() -> Point {
        Point::from_x(ModP::random(&mut thread_rng()))
    }

    pub fn is_zero(&self) -> Choice {
        self.z.is_zero()
    }

    pub fn normalize(&self) -> Point {
        let (z_inv, is_nonzero) = self.z.invert();
        if bool::from(is_nonzero) {
            Point {
                x: self.x * z_inv,
                z: ModP::ONE,
            }
        } else {
            Point {
                x: ModP::ONE,
                z: ModP::ZERO,
            }
        }
    }

    pub fn nonzero_to_bytes(&self) -> Option<[u8; 64]> {
        if bool::from(self.is_zero()) {
            return None;
        }
        Some(self.x.retrieve().to_be_bytes())
    }
}

impl fmt::Display for Point {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let p = self.normalize();
        write!(f, "{} : {}", p.x.retrieve(), p.z.retrieve())
    }
}

#[derive(Debug, PartialEq, Clone)]
pub struct MontgomeryCurve {
    pub a: Point, // projective coefficient, A = a.x/a.z
}

impl MontgomeryCurve {
    pub const fn new(a: ModP) -> MontgomeryCurve {
        MontgomeryCurve {
            a: Point::from_x(a),
        }
    }

    const fn projective(ax: ModP, az: ModP) -> MontgomeryCurve {
        MontgomeryCurve {
            a: Point { x: ax, z: az },
        }
    }

    pub fn twist(&self) -> MontgomeryCurve {
        MontgomeryCurve::new(self.a.x.neg())
    }

    pub fn normalize(&self) -> MontgomeryCurve {
        let a = self.a.normalize();
        MontgomeryCurve { a }
    }

    // Whether p corresponds to any rational point on the curve
    pub fn on_curve(&self, p: &Point) -> bool {
        let sq = p.x.square();
        let x3 = sq * p.x; // cube of p.x
        let x2 = sq * self.a.x;
        let rhs = x3 + x2 + p.x;
        // If rhs is square, point is only rational on this curve
        // If rhs is zero, point is rational on both curve and twist
        // If rhs is non-square, point is only rational on the twist
        //TODO: CT?
        let l = Integer::from_digits(
            &rhs.as_montgomery().as_limbs().map(|x| x.0.to_le()),
            Order::LsfLe,
        )
        .legendre(&P_GMP);
        // let l = jacobi_vartime(&rhs);
        l != -1
    }

    // Given points P, Q, and either P+Q or P-Q, computes P-Q or P+Q
    // P, Q must be distinct points, neither of which are the origin or infinity (zero)
    // Algorithm 1 in Costello & Smith
    fn add3(&self, p: &Point, q: &Point, pq: &Point) -> Point {
        let v1 = (p.x + p.z) * (q.x - q.z);
        let v2 = (p.x - p.z) * (q.x + q.z);
        let v3 = (v1 + v2).square();
        let v4 = (v1 - v2).square();
        Point {
            x: pq.z * v3,
            z: pq.x * v4,
        }
    }

    // Given P not equal to origin or infinity, computes P+P
    // Projective version of Algorithm 2 in Costello & Smith,
    // obtained by multiplying through by 4a.z
    fn double(&self, p: &Point) -> Point {
        // V1 in Costello & Smith
        let vplus = (p.x + p.z).square(); // (x + z)^2

        // V2 in Costello & Smith
        let vminus = (p.x - p.z).square(); // (x - z)^2

        let vdelta = vplus - vminus; // V1 in Costello & Smith
        let va = vminus * self.a.z;
        let va = va + va + va + va;

        let x = vplus * va; // vminus * vplus * denominator

        let vb = (self.a.z + self.a.z + self.a.x) * vdelta + va;

        Point { x, z: vb * vdelta }
    }

    // Given P not equal to origin or infinity, computes x([k]P)
    // Algorithm 4 in Costello & Smith, Montgomery ladder
    pub fn mult(&self, p: &Point, k: &ModP) -> Point {
        let k = k.retrieve();
        if bool::from(k.is_zero()) {
            return Point::zero();
        } else if k == U512::ONE {
            return p.clone();
        }
        let ell = k.bits() - 2;
        let mut x0 = p.clone();
        let mut x1 = self.double(p);

        // standard double-and-add, but adapted to x-only arithmetic
        for i in (0..=ell).rev() {
            let sum = self.add3(&x0, &x1, p);
            //TODO: CT?
            if bool::from(k.bit(i)) {
                x0 = sum;
                x1 = self.double(&x1);
            } else {
                x1 = sum;
                x0 = self.double(&x0);
            }
        }

        Point { x: x0.x, z: x0.z }
    }

    // Verify curve is nonsingular, i.e. A^2 != 4
    pub fn is_nonsingular(&self) -> bool {
        let a = self.a.normalize();
        const TWO: ModP = ModP::new(&U512::from_u8(2u8));
        if a.x == TWO {
            return false;
        } else if a.x + TWO == ModP::ZERO {
            return false;
        }

        return true;
    }

    // Algorithm 3 in Castryck et al.
    // Recursive subdivision to share computation of [(p+1)/ell_i]
    //   between different ell_i. lower is inclusive, upper is not
    fn is_supersingular_inner(
        &self,
        lower_idx: usize,
        upper_idx: usize,
        q: &Point,
        order_divisor: &mut ModP,
        z: &ModP,
    ) -> Option<bool> {
        if upper_idx - lower_idx > 1 {
            let midpoint = lower_idx + (upper_idx - lower_idx) / 2;

            // product of first half, second half of range of primes
            let lower_product: ModP = PRIMES[lower_idx..midpoint].iter().fold(ModP::ONE, |prod, p| prod * p);
            let upper_product: ModP = PRIMES[midpoint..upper_idx].iter().fold(ModP::ONE, |prod, p| prod * p);

            let qa = self.mult(q, &lower_product);
            let qb = self.mult(q, &upper_product);

            let upper = self.is_supersingular_inner(
                midpoint,
                upper_idx,
                &qa,
                order_divisor,
                &(lower_product * z),
            );
            if upper.is_some() {
                return upper;
            }

            let lower = self.is_supersingular_inner(
                lower_idx,
                midpoint,
                &qb,
                order_divisor,
                &(upper_product * z),
            );
            return lower;
        } else {
            let ell = PRIMES[lower_idx];
            let r = self.mult(&q, &ell); // ell_i * q = [p+1] * k
            if !bool::from(r.is_zero()) {
                return Some(false); // ell_i divides p+1 but not group order, so cannot be supersingular
            }

            // p+1 is multiple of group order but (p+1)/ell_i is not => ell_i divides group order
            if !bool::from(q.is_zero()) {
                *order_divisor = (*order_divisor) * ell;
            }

            // once the Hasse interval [(p+1) - 2*sqrt(p), (p+1) + 2*sqrt(2)] contains only one multiple
            // of the order_divisor, and order_divisor divides p+1, we can be sure the group order is p+1
            if order_divisor.retrieve() > *HASSE_INTERVAL {
                return Some(true);
            } else {
                return None;
            }
        }
    }

    // Algorithm 1 in Castryck et al.
    // Only valid for nonsingular curves
    pub fn is_supersingular(&self) -> bool {
        loop {
            let p = Point::random();
            const FOUR: ModP = ModP::new(&U512::from_u8(4u8));
            let p = self.mult(&p, &FOUR); // remove even factor from order
            if bool::from(p.is_zero()) {
                continue;
            }

            let mut order_divisor = ModP::ONE;

            let supersingular = self.is_supersingular_inner(
                0,
                PRIMES.len(),
                &p,
                &mut order_divisor,
                &FOUR,
            );
            if supersingular.is_some() {
                return supersingular.unwrap();
            }

            // failed to confirm or refute supersingularity
            // try again with another random point
        }
    }

    // Computes isogeny phi with kernel generated by point k, with order ell
    // Returns codomain curve of isogeny and image of point p under the isogeny
    // Algorithm from Castryck et al.
    pub fn isogeny(&self, k: &Point, ell: usize, p: &Point) -> (Point, MontgomeryCurve) {
        assert!(ell > 2);
        assert_eq!(ell % 2, 1);

        // inititalize computation of q = phi(p)
        let mut q = Point {
            x: ModP::ONE,
            z: ModP::ONE,
        };

        // shared factors for coefficient c computation
        let psub: ModP = p.x - p.z;
        let padd: ModP = p.x + p.z;

        // last three multiples of kernel generator k, kernel_points[i % 3] = [i+1]k
        let mut kernel_points: [Point; 3] = [k.clone(), self.double(k), Point::zero()];
        // let (x_i : z_i) = [i]k, and expand the product (z_i*w + x_i) as poly of w
        // c = [c_0, c_1, c_(ell-2), c_(ell-1)] are first two, last two coefficients of this expansion
        let mut c: [ModP; 4] = [ModP::ONE, ModP::ZERO, ModP::ZERO, ModP::ONE];

        // Since [i]P = [ell-i]P, we only compute up to ell/2 and then square the products
        for i in 0..(ell / 2) {
            if i > 1 {
                // [i]k = [i-1]k + k
                kernel_points[i % 3] =
                    self.add3(&kernel_points[(i - 1) % 3], k, &kernel_points[(i - 2) % 3]);
            }

            // Xi/Zi = [i+1]P
            let xi: &ModP = &kernel_points[i % 3].x;
            let zi: &ModP = &kernel_points[i % 3].z;

            // multiply polynomial by (Zi*w + Xi)
            let a = c[0] * zi;
            let b = c[1] * xi;
            c[1] = a + b;
            c[0] *= xi;

            let a = c[2] * zi;
            let b = c[3] * xi;
            c[2] = a + b;
            c[3] *= zi;

            // Inner portion of equation 17 in Costello & Hisil, computes:
            // a = (X-Z)*(Xi+Zi), b = (X+Z)*(Xi-Zi)
            // Q.x *= (a + b)
            // Q.z *= (a - b)
            let ikadd = xi + zi;
            let iksub = xi - zi;
            let a = psub * ikadd;
            let b = padd * iksub;
            let s = a + b;
            let t = a - b;

            q.x *= s;
            q.z *= t;
        }

        // Compute Q = phi(P) via Equation 17 from Costello & Hisil
        q.x = p.x * q.x.square();
        q.z = p.z * q.z.square();

        // square the polynomial
        c[1] *= c[0];
        c[1] += c[1];
        c[0] = c[0].square();

        c[2] *= c[3];
        c[2] += c[2];
        c[3] = c[3].square();

        // Codomain coefficient formula from Castryck et al.
        // A'.x = A.x*c[0]*c[ell-1] - A.z*3*(c[0]*c[ell-2] - c[1]*c[ell-1])
        // A'.z = A.z*c[ell-1]^2
        let a = self.a.x * c[0] * c[3];
        let b = (c[0] * c[2]) - (c[1] * c[3]);
        let b = self.a.z * b;
        let b3 = b + b + b;

        let ax = a - b3;
        let az = self.a.z * c[3].square();
        let codomain = MontgomeryCurve::projective(ax, az);
        (q, codomain)
    }

    pub fn to_bytes(&self) -> Option<[u8; 64]> {
        self.a.nonzero_to_bytes()
    }

    pub fn from_be_slice(a: &[u8]) -> MontgomeryCurve {
        MontgomeryCurve::new(ModP::new(&U512::from_be_slice(a)))
    }
}

#[cfg(test)]
mod tests {
    use crate::constants::BASE_CURVE;

    use super::*;

    #[test]
    fn on_curve() {
        let a = ModP::new(&U512::from_be_hex("57164DAD2DAA6B17538CC28E418D0B93540024EC5F7038951049142A2FA46F030D7E5247B792A894FF526D7126DCB9CDEF42704493B6F8109CC5B127FD6F4888"));
        let e = MontgomeryCurve::new(a);

        let p = Point::from_x(ModP::new(&U512::from_be_hex("39F65EDE480BF8C5E5ACA9CE8EAC7EC98B8C02E3768B83444D77F06961B50B9CDCD8A3644D244624766ABECED59F881BA06B04033B6FA7652396BA8798A16CEA")));
        assert!(e.on_curve(&p));

        let p = Point::from_x(ModP::new(&U512::from_be_hex("4ACCB8B58B35F4D9787314CB062D264A5CD43EC672B48CEAD6FE63FD49A94CC36912F751EAE0D262F1584DA663F2A3A18506EF9F4B444D7F40AAC0A7E5838869")));
        assert!(e.on_curve(&p));

        let p = Point::from_x(ModP::new(&U512::from_be_hex("5727FDAA0E1070CFF054606C531F4BE7B6D55B9B2CC2F343C19306C76BFA01D247BBDC9F05E5AFD1D08453C532F1E733CC3419C868167AF4F03AC860A90E258F")));
        assert!(!e.on_curve(&p));

        let p = Point::from_x(ModP::new(&U512::from_be_hex("63576655EC3AA520BF7A2B022D5253C18E8676C03EF81FA05030B2A4509F2E4C11A37DDDB06606A6CA94DAD30DD0876E300FD9AEAEF2B075DD5CE6AD74255B13")));
        assert!(!e.on_curve(&p));
    }

    #[test]
    fn add3() {
        let m = ModP::new(&U512::from_be_hex("54C8A0DADC3C7B204FDF48616AF757968326DC25866F018424FCD27D45C809CAC0D3F58553D6CB42704819843C67406977C51CD790BE78350FADAB6CB72AFA8D"));
        let mut e = MontgomeryCurve::new(ModP::ZERO);
        e.a.x *= m;
        e.a.z *= m;

        let p = Point::from_x(ModP::new(&U512::from_be_hex("2C269236F8C147E177F0A3B5C95EA91A4C04DE1BBA1FFC84622CE367805C8551A7D628FFEA33900BA7F6F90F65AA0EAFEA189CEC225A1730B7F1E28FF2C253B4")));
        let q = Point::from_x(ModP::new(&U512::from_be_hex("258DD34D9A793CE8396185303C65012683D12744E672B2195B38D7A2DF26242024C3E8740CCAFD3071388C488266E2BBCBA9F5B3F2D4252CBCDB4EF69DD58EDB")));
        let pq_sum = Point::from_x(ModP::new(&U512::from_be_hex("55DA2197E667E35145692AB34A06D9C62A379533C8430B98295F3F0A6828682BC7046F3DF84760F6C22CBBEA55B67EF83C53AB1CF103DE0CF0774719D9038C46")));
        let pq_diff = Point::from_x(ModP::new(&U512::from_be_hex("5A54CDD887887905BB57A2FB2C5D108A361052CAF9B80E93812C69E5BF5EA6C9869F43748FB12937C5A91C999AAF4A5D0CBA6A0B9B67A8339F15DA786626F7F3")));

        assert_eq!(e.add3(&p, &q, &pq_sum).normalize(), pq_diff);
        assert_eq!(e.add3(&p, &q, &pq_diff).normalize(), pq_sum);

        let a = ModP::new(&U512::from_be_hex("47D112C8D0BBF39D1983F677BE0CD423445C8BACA91B516EB3350F1CB95FFB454F4B0C18CE2EA540CE7B0932B951B365511CDBB82458DCA4D0ABBA04DB00D84D"));
        let m = ModP::new(&U512::from_be_hex("52C521933A01AD67352ABAEE2BB6FDB4025BA653A1B6C5C8B939B5647EF56A8111640A7717FEEB38967FE1F7B653384D1E4C5E4DAA76686F97CA4DC9E405C543"));
        let mut e = MontgomeryCurve::new(a);
        e.a.x *= m;
        e.a.z *= m;

        let p = Point::from_x(ModP::new(&U512::from_be_hex("55D5F791250BF6F55EF6DCDE2F1249DB59DAEA4FDAB940A53FC796293B137D58FC0888CF466C523C21EED52CDA1554705E3ECE421FDB0653E052DE3009F60C29")));
        let q = Point::from_x(ModP::new(&U512::from_be_hex("38AC0C3673F3DE0E5C049EA68C97CE5261EEE5EC6BB02F568E0D3B70A40E139F22C64995A40B1B5F90D3AFF943557C0E14BB8227C477FC7E822C5787450198DB")));
        let pq_sum = Point::from_x(ModP::new(&U512::from_be_hex("5C58A818F02D09EF5DC1CAABEC133D07D3D9EFB98984B1DD9B7636D4A4FBC6B47FBF51FCDB71570CE00C723F4083CF4E70F1467077FC99B0AB746004ECC445BD")));
        let pq_diff = Point::from_x(ModP::new(&U512::from_be_hex("174B7233B392E09C38F7DB5DA5D352585FBE4D16CBD626404D894CA3407AB399400FD9FFE0FDF1C76ADB7F11955F88F45805EBDD727271F5465BEB39D96E6782")));

        assert_eq!(e.add3(&p, &q, &pq_sum).normalize(), pq_diff);
        assert_eq!(e.add3(&p, &q, &pq_diff).normalize(), pq_sum);
    }

    #[test]
    fn double() {
        let m = ModP::new(&U512::from_be_hex("09F4B54FF7BCF319672BEDF8D865F241C27D52BE5D70F8E4BD18806AF9E5BFF3FC85635574EC6D8513679E8F1BBE290672AC5AB4A0106D05C10DD74B758F2589"));
        let mut e = MontgomeryCurve::new(ModP::ZERO);
        e.a.x *= m;
        e.a.z *= m;

        let p = Point::from_x(ModP::new(&U512::from_be_hex("24FECD3460B2F53C30EC39BE2E4CE435EAC8EE1E5E57CC2A022D9C885F941346D78DFA27AF0BE7CBFE913E461310BB5713A85D7AC02AAD5B8E802D2F7D045199")));
        let q = Point::from_x(ModP::new(&U512::from_be_hex("50809AE8A1F22F00D2D655E4E6E7580A1C868E59AE4BE4E98A7AF3FD3BB0A26C24724B6E0C98EFA529C914E7F5DFA89D1CBBCEF8E2DB859AB9538A7787B10EFB")));
        let two_p = Point::from_x(ModP::new(&U512::from_be_hex("08D21FCB978BE16F2A6A9A5753418098BE7CE928B417E4672B448F99029BF4F9D771C6AE5356889A0A7F0DD30AC467343CEDFD25BDE221C8B82B9D208CB3F2C1")));
        let two_q = Point::from_x(ModP::new(&U512::from_be_hex("02C862DCE53D5BE8C3FDDBC31152CAF717D52579033EF877E8DCAD1BF333C3B590D045023BD3B365D0FB0EDDAFB2D74871BA618EAB37C285B91CD83CF395F848")));

        assert_eq!(e.double(&p).normalize(), two_p);
        assert_eq!(e.double(&q).normalize(), two_q);

        let a = ModP::new(&U512::from_be_hex("20A24B2BE41B432B63B7C485E95F781E344744881BB3C9C36A10995DFB62592FFFDF6EAEA75F632894CC2D30F0C273AF4002D580E272FC2B4EB2C46730A6D7CD"));
        let m = ModP::new(&U512::from_be_hex("3DB6207324441ECC3884C3865E972993CA0C239CAABDC1E4BFEB70E73CA776A83C8AF9E383B93CEDE3E6BA05503D422EEE7FAA5484F4711E69939FA3FDF7DD97"));
        let mut e = MontgomeryCurve::new(a);
        e.a.x *= m;
        e.a.z *= m;

        let p = Point::from_x(ModP::new(&U512::from_be_hex("228B3C1B541F1C50E6E46545E9E3E26F747EC6D18D436CD02DA23B8AF7A37FD81F1132833AC5F96EBD3AF861F2CD83EAF429017B859EB55D6A7BAF514849B1BB")));
        let q = Point::from_x(ModP::new(&U512::from_be_hex("3F53508C3824E5B18630C63E35FD312E39B36345E7BA855D1AE5DA044682327A8FBE7627C46E707CD4E5F87BD0EB191D366DC0A7D7F2ECDFC7E6F9ECFC6A55D2")));
        let two_p = Point::from_x(ModP::new(&U512::from_be_hex("29306C62A1CFDDDEFE22DB90F5C3BF81DDA85BF6F25B45F37D6ECA62D633066CE47D707B9698EED38F1C2FD9D60E13EF33B712603FEE555A2EEFF68177B50A47")));
        let two_q = Point::from_x(ModP::new(&U512::from_be_hex("0E417B11BA938C20C957E22F787757DD7B138C633DF90F529DD0EC490A7FA49E847644AE5728E981DDCFC24D37B93EF169502C66FDBDE16D6D0A986DB0640738")));

        assert_eq!(e.double(&p).normalize(), two_p);
        assert_eq!(e.double(&q).normalize(), two_q);
    }

    #[test]
    fn mult() {
        let a = ModP::new(&U512::from_be_hex("0F65404A8762A1310462B75091BC880B90F5CC563F59FBC0895112E35F27E6DE9FBC44F95037C11D1C380BB6FBA089A3ABE20A6B2A89332BAF14CE916900268E"));
        let m = ModP::new(&U512::from_be_hex("2CC1305DAF90CF3849C4AE87B8C5F969ABDA84BDCF9356019318433CD87AA3F0BAC9455F2E2E6B233D78A6ADB31C2D233AE4818F34C837A4BE7B63570180714C"));
        let mut e = MontgomeryCurve::new(a);
        e.a.x *= m;
        e.a.z *= m;

        let p = Point::from_x(ModP::new(&U512::from_be_hex("1D5FEE8C54EAE6791B020BDE73090A0CDE5AAF6472F6C251C3BB46EA2DADC6E412CFE135AA6B65E11D1AB4B9F612918D2C050B6C105708BC2814A41BB71A1314")));
        let q = Point::from_x(ModP::new(&U512::from_be_hex("638D2CC4889113F79C022F08BC6EAD9A95EBD83CFD914BE639E2E6CCDB70FC6876CB9DFF05AC06698FE172B327C59FA605645B47D07BB15EFA2F7247A03EE84B")));

        let k = ModP::ZERO;
        assert_eq!(e.mult(&p, &k).normalize(), Point::zero());
        assert_eq!(e.mult(&q, &k).normalize(), Point::zero());

        let k = ModP::ONE;
        assert_eq!(e.mult(&p, &k).normalize(), p);
        assert_eq!(e.mult(&q, &k).normalize(), q);

        let k = ModP::new(&U512::from_be_hex("1AFD4893406045345813F35A4E8519DF1D126D6B0FA04259312890AC1011592DC33B75EC8F3427C5BF620F095A06E0E28557389CA1AE41DC6D671D256E132924"));
        let kp = Point::from_x(ModP::new(&U512::from_be_hex("36C63051A245BEAF9C715197C7EB37C9F199ECB72C5ADE3376E0290E8E95E2AF477B8CDC2588DB041B8FD455B62FB2B26CC066B6EE85072F0BA1F57CCE119979")));
        let kq = Point::from_x(ModP::new(&U512::from_be_hex("4105472DBE22538577C4144B0DFBA98977F74DF8EC235E5876484803ED8C4C34CB72534E834AB72D86138A1F5B01A006DF1AAA6030834403B811D2F6FF74BEB7")));
        assert_eq!(e.mult(&p, &k).normalize(), kp);
        assert_eq!(e.mult(&q, &k).normalize(), kq);

        let k = ModP::new(&U512::from_be_hex("0D4B1AEE28AB783B589255E0C0B6983DE8F5B3ACEE733C36B9A8B7B362A3AE94B0BB71595D324F8065A3B32A81E3DCD5FC26CB02CA763ABF7CCAD55F96DB11EC"));
        let kp = Point::from_x(ModP::new(&U512::from_be_hex("145730A5A88421AFFD62098BAE43D698E269F011D926BB37060F2699D367CDD0213664632B963EEEB649F31B2394C3E0B1264EC021CD62215234F1A6DACF4634")));
        let kq = Point::from_x(ModP::new(&U512::from_be_hex("3F18A2D2C8814FD94B1AF5888073F4BE0EDBF19E11A86F953FFC8C20A028CBA6E6E7201CB9277A950EEB8543D05D758B6CC1F425C4BAB16963E813D095121E49")));
        assert_eq!(e.mult(&p, &k).normalize(), kp);
        assert_eq!(e.mult(&q, &k).normalize(), kq);
    }

    #[test]
    fn is_nonsingular() {
        let m = ModP::new(&U512::from_be_hex("1c50b5fd8d624c6fc1b4cf7df10d56bcd193fa0a32e1097946639b8ea84a819f4506d3ff77fe9c723cb26f0ba94e5775a8ea103a8ea777e241b4c6b5d9a7503a"));
        let mut e = BASE_CURVE;
        e.a.x = (&e.a.x) * m;
        e.a.z = (&e.a.z) * m;
        assert!(e.is_nonsingular());

        let m = ModP::new(&U512::from_be_hex("4353650d3a92105abfdd02136d6ff90a8c1833a9ecb5d048fee7b444e53d5d59cde0e4829b767a9014fa22cef7df61d6d9e2b042a100f9f053b50fa1b4d569ad"));
        let mut e = MontgomeryCurve::new(ModP::new(&U512::from_u8(2u8)));
        e.a.x = (&e.a.x) * m;
        e.a.z = (&e.a.z) * m;
        assert!(!e.is_nonsingular());

        let minus_two: ModP = -ModP::new(&U512::from_u8(2u8));
        let m = ModP::new(&U512::from_be_hex("5d4d97e81b06c5344156f9da5e885e03fec07b0bad7cb707efa0c56ad22c528705b1efd093358ceb6f8e30938b76eb55b0e99c249d2f842810c89677d0d12c1e"));
        let mut e = MontgomeryCurve::new(minus_two);
        e.a.x = (&e.a.x) * m;
        e.a.z = (&e.a.z) * m;
        assert!(!e.is_nonsingular());

        let a = ModP::new(&U512::from_be_hex("23a1e90cebfe1fafa16a65820c2f19818d031c1832d5905349f24192ec5303eeffaebfd4652bea5c845a8353d00d2b084cad074636f1d8c4aabd4907e9e3a5af"));
        let m = ModP::new(&U512::from_be_hex("434a044d8e844abd947d9a33835e32a6bfb1056a6a6823821239dc45951fcd9903d8bf5aaa19ee4d41ee8c2f0a7c032b0ef15ac2cec27543c71b856333075515"));
        let mut e = MontgomeryCurve::new(a);
        e.a.x = (&e.a.x) * m;
        e.a.z = (&e.a.z) * m;

        assert!(e.is_nonsingular());

        let a = ModP::new(&U512::from_be_hex("20a24b2be41b432b63b7c485e95f781e344744881bb3c9c36a10995dfb62592fffdf6eaea75f632894cc2d30f0c273af4002d580e272fc2b4eb2c46730a6d7cd"));
        let m = ModP::new(&U512::from_be_hex("29adb021a7afc685859661441c4f1854f1b1a388a9bed73c0dbc5171dc1f3a88c5147b4d637f34ccc4a90aa2379723df024ae36f160391e220fe23e9c11f001c"));
        let mut e = MontgomeryCurve::new(a);
        e.a.x = (&e.a.x) * m;
        e.a.z = (&e.a.z) * m;

        assert!(e.is_nonsingular());
    }

    #[test]
    fn is_supersingular() {
        // Supersingular curves
        let m = ModP::new(&U512::from_be_hex("57164dad2daa6b17538cc28e418d0b93540024ec5f7038951049142a2fa46f030d7e5247b792a894ff526d7126dcb9cdef42704493b6f8109cc5b127fd6f4888"));
        let mut e = BASE_CURVE;
        e.a.x = (&e.a.x) * m;
        e.a.z = (&e.a.z) * m;
        assert!(e.is_supersingular());

        let a = ModP::new(&U512::from_be_hex("47d112c8d0bbf39d1983f677be0cd423445c8baca91b516eb3350f1cb95ffb454f4b0c18ce2ea540ce7b0932b951b365511cdbb82458dca4d0abba04db00d84d"));
        let m = ModP::new(&U512::from_be_hex("5611c008258a404f2bca6eb568ccde8993fd983e3bfb99415eeb8bdca00640a42d6184646ecb8f78be29bb2142f3481eb70a031395d663e2d4dfa222f95ac19b"));
        let mut e = MontgomeryCurve::new(a);
        e.a.x = (&e.a.x) * m;
        e.a.z = (&e.a.z) * m;
        assert!(e.is_supersingular());

        // let a = BigUint::from_str("806328403495992267109149450868966636991037459819762153296867249317594469447447955344137556739687998616446227640966045492819498968091824989742441223431822"));
        // let m = BigUint::from_str("682589759858867006019945921981683911930844413103153096751081265553732323158453661898923474704061121351403714997566638936235758098864420090141210600541605"));
        // let mut e = MontgomeryCurve::new(a);
        // e.a.x = (&e.a.x) * m;
        // e.a.z = (&e.a.z) * m;
        // assert!(e.is_supersingular());

        // let a = BigUint::from_str("1709179145736732726190467057772287838888749898111981175668039410052353567995718751503892000081725932287177755925719395202935581874824412236400752119437261"));
        // let m = BigUint::from_str("1372299461168028207574437529954930118298843176446652665420970385704081105088164666433577075908289791477766172637812347564025927002291465951964304934275749"));
        // let mut e = MontgomeryCurve::new(a);
        // e.a.x = (&e.a.x) * m;
        // e.a.z = (&e.a.z) * m;
        // assert!(e.is_supersingular());

        // Ordinary curves
        let a = ModP::new(&U512::from_be_hex("250a0b2656eccdc6d884a9e781edbc782ab90290432c33ddb339d357358ffff6b83b0facbf99e1c1ab1b437368e95082593c8528c8d24871f4f1bdb88ad11660"));
        let m = ModP::new(&U512::from_be_hex("4d3265cdc86303738bb64ae6bce527757d19bbec94dd57933be1654055e148d29bdb4b09c0ad90f23d423e3b5d22673486edce46a857a748d0f6211cdf17c0d2"));
        let mut e = MontgomeryCurve::new(a);
        e.a.x = (&e.a.x) * m;
        e.a.z = (&e.a.z) * m;
        assert!(!e.is_supersingular());

        // let a = BigUint::from_str("1866223479810078799335388829853753029730709787935267021757606377075016878572548078343651839809404146780592040359635385767112537414085805658738687484339631"));
        // let m = BigUint::from_str("3640326944036693676475257727036086990099851523941743651295162127638831471205355580550155137235405452465381206319509761479787085502435824357185503728478626"));
        // let mut e = MontgomeryCurve::new(a);
        // e.a.x = (&e.a.x) * m;
        // e.a.z = (&e.a.z) * m;
        // assert!(!e.is_supersingular());

        // let a = BigUint::from_str("2661785568965525579235733734106882226803961945735665715347090419442451011470046364314754336833920717330106734294787860943873523408164640783182492955709999"));
        // let m = BigUint::from_str("683201296837654516155062671807374643518145895323214941072784881241198594642968581707421521295675835749259753224602137584501907664571532242696362741510003"));
        // let mut e = MontgomeryCurve::new(a);
        // e.a.x = (&e.a.x) * m;
        // e.a.z = (&e.a.z) * m;
        // assert!(!e.is_supersingular());
    }

    #[test]
    fn isogeny() {
        let a = ModP::ZERO;
        let b = ModP::new(&U512::from_be_hex("53BAA451F759835A01933C76BC58C0C203A9B6B02F7F086B30C3469A8452750AAECA8A4F7C26BFF43876F4510F405F4D2A006635D89A42D327D9A2E8C00BF340"));
        let order: usize = 3usize;

        let k = Point::from_x(ModP::new(&U512::from_be_hex("22B668C942BF7D5F5DF869A215F7E9463A0A873CFE2953721F129EC98B8123A8E62DF0D1F100AA92F4C6C8552AD62C42C11DB1AE8540F46ADC16D8939808553A")));
        let p = Point::from_x(ModP::new(&U512::from_be_hex("0A3A72458C434F22FD1F2B441C3BAD38C0C069872F69372A43E818126CFF49DC3CA63E87BC5F0443201F9DA03EFE8DA618C4D207954D40F774A923CBC11F2CA7")));
        let im_p = Point::from_x(ModP::new(&U512::from_be_hex("1ED168610F98DC95AAB55E2B067E92B32AF0A436A73EF7142F31BC3CBE2A532F8D51061DA110C5EB01FEC1838C6D0AA3B643D90181AAA3184CF02ABB20ECFB2A")));

        let ea = MontgomeryCurve::new(a);

        let (xp, xb) = ea.isogeny(&k, order, &p);
        let xb = MontgomeryCurve::new(xb.a.normalize().x);
        assert_eq!(xb.a.x, b);
        assert_eq!(im_p.normalize().x, xp.normalize().x);

        let a = ModP::new(&U512::from_be_hex("48211766D23E629D22C38ED44B3D8A02622B7022E5CE2CE5CCF7CDD4F901213AE61B00371E74AD24C9F71C59C0B0269287B36EC9652F4ACC421B8975C8C9EE4F"));
        let b = ModP::new(&U512::from_be_hex("20B68C844B20BBA2271497C8ECD471D2EB0E3640A3D238F142C13C3C86BDF9D2F758186586740B2A15F9709E18F93F7894704B23CCBB533AC8AD2F1031AE309B"));
        let order: usize = 5usize;
        let k = Point::from_x(ModP::new(&U512::from_be_hex("409548FAF7B5117391A5AD4D1202CA9EE096D69F44188441796F2ACED23C0C21DA29C9286AD5A46636CE1E41F9F54CEF4F453F7EFFCD595E168CC519DD68EA51")));
        let p = Point::from_x(ModP::new(&U512::from_be_hex("3B2DC0FD5EE8C65F43DAD597D8C48C32138A9FE4A1008802D5CED33523731EB432469E2D7F2276625E3DF38566576180E559E1C13D5F9565696A6D0D83830FF4")));
        let im_p = Point::from_x(ModP::new(&U512::from_be_hex("3BF01DE995EB675B0C2303367BC0FC3F82AC3D7123F842DEC8DE1E34F6FE14FBFDC1BDF203914BF7F6C52A7AF66CA745D4C682A4D1C9F40D1D5CFB066BB46B2D")));
        let ea = MontgomeryCurve::new(a);

        let (xp, xb) = ea.isogeny(&k, order, &p);
        let xb = MontgomeryCurve::new(xb.a.normalize().x);
        assert_eq!(xb.a.x, b);
        assert_eq!(im_p.normalize().x, xp.normalize().x);

        let a = ModP::new(&U512::from_be_hex("48912381B13014DF7E10F242424DFE6D43860ED48A2913843A45E75E15615849B2E2C8191E6CEF70A931E20883E8B59B87046926B8E534DCA88722A8E204496C"));
        let b = ModP::new(&U512::from_be_hex("4E7D723C463A2F779721CD1F53CB1F2F3F9ECEB60E1831A2DDC665687C1F7BD1B479670592F4967DFAD3F9675B6229ED2B4ECFC2AAE33258DAF6A1B5CCBF2B78"));
        let order: usize = 173usize;
        let k = Point::from_x(ModP::new(&U512::from_be_hex("54C8CDA4F5B40B4DD5EF9011AFA313195A68106114B157B53270CB1005C8338E4CE00C826ECEE406027F383FA1D5037DBB81D92E4203B4092B9C3D20A32D49A8")));
        let p = Point::from_x(ModP::new(&U512::from_be_hex("55E59A6BB770F1477F38E747D4C45F61CDA4D068736398DCB7C3A6B872208E6BA55FA42377A4B3EB25AEF4D0CE59C91A1D3A291B87700FFFE21805A7DEED199B")));
        let im_p = Point::from_x(ModP::new(&U512::from_be_hex("5BBC4E76109DAF1BCC7A597C78DCA56C1645CA6C72859B3F316F972054BF200C0F8059E2CCD9B1886F7518230CC5A75E210A3A5C07D843FF79BD832B675E1BD3")));
        let ea = MontgomeryCurve::new(a);

        let (xp, xb) = ea.isogeny(&k, order, &p);
        let xb = MontgomeryCurve::new(xb.a.normalize().x);
        assert_eq!(xb.a.x, b);
        assert_eq!(im_p.normalize().x, xp.normalize().x);

        let a = ModP::new(&U512::from_be_hex("4F80CF43DC32028D21AF9A4596C7067352C764156B62056D1DBC2A528E367DAC0BC65E453F01864DF53E0A775E064988EDD71034EBE1D5B95C1235F11DC6ACFD"));
        let b = ModP::new(&U512::from_be_hex("1AFAA394A786BAEB11895EE8A455AE6A6872C74C9D2F0F47773AAB2FD1481BACC7695E7B81177C643C054D3BD36268F5ABD7AE225EEBEFA531F9153F532E4BA5"));
        let order: usize = 83usize;
        let k = Point::from_x(ModP::new(&U512::from_be_hex("425A0D7407BF49078B071367E138506CDF3CF5C5231384524F9C62C7E84BF1536C47B5AC7B981BD1D8B8A4FF5ED0A75471F0D80ED1515CE18C2D31780929FD58")));
        let p = Point::from_x(ModP::new(&U512::from_be_hex("1C33BAADF7E34ACA1AFE98CEF02BA3948B0A09A2996BC9BC2C28A4E33A4943E8FA63370DF59C0CDA3C1D943473E50B2D4334DEE8263F6CF6450619460BF09DC4")));
        let im_p = Point::from_x(ModP::new(&U512::from_be_hex("46C02021FCC07A06A4C41958C3C40EB31DF7D10947C1021FE2638A9DADB7C8792D8EC0271FE63DFDE0BF6E1B4D44E550A9606DC91541DC15263292469892BD8D")));
        let ea = MontgomeryCurve::new(a);

        let (xp, xb) = ea.isogeny(&k, order, &p);
        let xb = MontgomeryCurve::new(xb.a.normalize().x);
        assert_eq!(xb.a.x, b);
        assert_eq!(im_p.normalize().x, xp.normalize().x);

        let a = ModP::new(&U512::from_be_hex("38C29C0734D6DD3278A7DCDB797A4D8D1E6C41ADAC0FA768EE3BF9EE8FBC3EECE077F09DE1644631226B822B3BA01868DE4A9E7603DADA2D024BB5B16E083020"));
        let b = ModP::new(&U512::from_be_hex("2FDC10B37A9570DFA25AEE9482802ECDC60E7D5D47B6E06AC5C3114BA70DDF38C6D820DFAD5A126794DC0CCF78A3BB91283DEAE9D6B540EB45506934D5C145B7"));
        let order: usize = 149usize;
        let k = Point::from_x(ModP::new(&U512::from_be_hex("5E9D1A6638D9610AE0568BC36A483E512E3AC582C45E79A5388D0C213F7315052B4B74784F468E9D5CE9EA882D9511AF4A7B92E7CDEC4D5AE22D32D8B9F805DE")));
        let p = Point::from_x(ModP::new(&U512::from_be_hex("35EE2441FB15C6C330837FA2950A9860C33A14E6847D78DE6EB62FF85291477CEB7E69CE825B88637283A87379AC17D3A1E319A2D95172CFFC31FE6380C54749")));
        let im_p = Point::from_x(ModP::new(&U512::from_be_hex("3A12C8F02C85892393291F5860DB7B8C86C198FE89B44B165A91E9C05F185896C036B64331A418347706C6D124B73AECE248925112F207E3E53114FEECE14545")));
        let ea = MontgomeryCurve::new(a);

        let (xp, xb) = ea.isogeny(&k, order, &p);
        let xb = MontgomeryCurve::new(xb.a.normalize().x);
        assert_eq!(xb.a.x, b);
        assert_eq!(im_p.normalize().x, xp.normalize().x);
    }
}
