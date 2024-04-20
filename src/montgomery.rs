use num_bigint::{BigUint, RandBigInt};
use num_modular::{ModularCoreOps, ModularPow};
use num_traits::{One,Zero};
use std::fmt;

use crate::{classgroup::{HASSE_INTERVAL, P, PRIMES}, modular};

#[derive(Clone, Debug, PartialEq)]
pub struct Point {
    pub x: BigUint,
    pub z: BigUint,
}

impl Point {
    pub fn from_x(x: BigUint) -> Point {
        return Point {
            x, z: One::one()
        };
    }

    pub fn zero() -> Point {
        return Point {
            x: BigUint::one(),
            z: BigUint::zero(),
        };
    }

    pub fn is_zero(&self) -> bool {
        return self.z.is_zero();
    }

    pub fn normalize(&self) -> Point {
        if self.z.is_zero() {
            return Point{x: BigUint::one(), z: BigUint::zero()};
        } else {
            let z_inv = modular::inverse(self.z.clone(), &P);
            let x = (&self.x).mulm(&z_inv, &P);
            return Point{x, z: BigUint::one()};
        }
    }
}

impl fmt::Display for Point {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let p = self.normalize();
        return write!(f, "{} : {}", p.x, p.z);
    }
}

// Elliptic curve in Montgomery form y^2 = x^3 + ax^2 + x
//     See Costello & Smith, https://eprint.iacr.org/2017/212.pdf for details
//     Also see CSIDH paper, Castryck et al., https://csidh.isogeny.org/csidh-20181118.pdf
// Curves are defined over the finite field with P elements
// TODO: implement Montgomery arithmetic for BigUint?
pub struct MontgomeryCurve {
    pub a: Point, // projective coefficient, A = a.x/a.z
}

impl MontgomeryCurve {
    pub fn new(a: BigUint) -> MontgomeryCurve {
        return MontgomeryCurve { a: Point::from_x(a) };
    }

    // Given points P, Q, and either P+Q or P-Q, computes P-Q or P+Q
    // P, Q must be distinct points, neither of which are the origin or infinity (zero)
    // Algorithm 1 in Costello & Smith
    fn add3(&self, p: &Point, q: &Point, pq: &Point) -> Point {
        let square = BigUint::from(2u32);

        let v0 = (&p.x).addm(&p.z, &P);
        let v1 = (&q.x).subm(&q.z, &P);
        let v1 = &v1.mulm(&v0, &P);

        let v0 = (&p.x).subm(&p.z, &P);
        let v2 = (&q.x).addm(&q.z, &P);
        let v2 = v2.mulm(&v0, &P);

        let v3 = &v1.addm(&v2, &P);
        let v3 = &v3.powm(&square, &P);

        let v4 =  &v1.subm(&v2, &P);
        let v4 = &v4.powm(&square, &P);

        return Point {
            x: (&pq.z).mulm(v3, &P),
            z: (&pq.x).mulm(v4, &P),
        }
    }

    // Given P not equal to origin or infinity, computes P+P
    // Projective version of Algorithm 2 in Costello & Smith,
    // obtained by multiplying through by 4a.z
    fn double(&self, p: &Point) -> Point {
        let square = BigUint::from(2u32);

        let v1 = (&p.x).addm(&p.z, &P);
        let v1 = v1.powm(&square, &P);

        let v2 = (&p.x).subm(&p.z, &P);
        let v2 = v2.powm(&square, &P);

        let v3 = (&v1).subm(&v2, &P);

        let v2 = (&v2).mulm(&self.a.z, &P);
        let v2 = (&v2).addm(&v2, &P);
        let v2 = (&v2).addm(&v2, &P);

        let x = (&v1).mulm(&v2, &P);

        let v4 = (&self.a.z).addm(&self.a.z, &P);
        let v4 = (&self.a.x).addm(&v4, &P);
        let v4 = (&v4).mulm(&v3, &P);

        let v4 = (&v4).addm(&v2, &P);

        return Point {
            x: x.clone(),
            z: (&v4).mulm(&v3, &P),
        }
    }

    // Given P not equal to origin or infinity, computes x([k]P)
    // Algorithm 4 in Costello & Smith, Montgomery ladder
    pub fn mult(&self, p: &Point, k: &BigUint) -> Point {
        if k.is_zero() {
            return Point::zero();
        } else if k.is_one() {
            return p.clone();
        }

        let ell = k.bits();
        let mut x0 = (*p).clone();
        let mut x1 = self.double(p);

        assert!(ell >= 2);

        // standard double-and-add, but adapted to x-only arithmetic
        for i in (0..=(ell-2)).rev() {
            let sum = self.add3(&x0, &x1, &p);

            if k.bit(i) {
                x0 = sum;
                x1 = self.double(&x1);
            } else {
                x0 = self.double(&x0);
                x1 = sum;
            }
        }

        return Point {
            x: x0.x,
            z: x0.z
        }
    }

    // Verify curve is nonsingular, i.e. A^2 != 4
    pub fn is_nonsingular(&self) -> bool {
        let a = self.a.normalize();
        if a.x == BigUint::from(2u32) {
            return false;
        } else if (&a.x).addm(BigUint::from(2u32), &P) == BigUint::zero() {
            return false
        }

        return true;
    }

    // Algorithm 3 in Castryck et al.
    // Recursive subdivision to share computation of [(p+1)/ell_i]
    //   between different ell_i. lower is inclusive, upper is not
    fn is_supersingular_inner(&self, lower_idx: usize, upper_idx: usize, q: &Point, order_divisor: &mut BigUint, z: &BigUint) -> Option<bool> {
        if upper_idx - lower_idx > 1 {
            let midpoint = lower_idx + (upper_idx - lower_idx)/2;
            let lower_product: BigUint = PRIMES[lower_idx..midpoint].iter().product(); // product of first half
            let upper_product: BigUint = PRIMES[midpoint..upper_idx].iter().product(); // product of second half

            let qa = self.mult(q, &lower_product);
            let qb = self.mult(q, &upper_product);

            let upper = self.is_supersingular_inner(midpoint, upper_idx, &qa, order_divisor, &lower_product.mulm(z, &P));
            if upper.is_some() {
                return upper;
            }

            let lower = self.is_supersingular_inner(lower_idx, midpoint, &qb, order_divisor, &upper_product.mulm(z, &P));
            return lower;
        } else {
            let ell = BigUint::from(PRIMES[lower_idx]);
            let r = self.mult(&q, &ell); // ell_i * q = [p+1] * k
            if !r.is_zero() {
                return Some(false); // ell_i divides p+1 but not group order, so cannot be supersingular
            }

            // p+1 is multiple of group order but (p+1)/ell_i is not => ell_i divides group order
            if !q.is_zero() {
                let z: &BigUint = order_divisor;
                *order_divisor = z.mulm(&ell, &P);
            }

            // once the Hasse interval [(p+1) - 2*sqrt(p), (p+1) + 2*sqrt(2)] contains only one multiple
            // of the order_divisor, and order_divisor divides p+1, we can be sure the group order is p+1
            if *order_divisor > *HASSE_INTERVAL {
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
            let x = rand::thread_rng().gen_biguint_below(&P);
            let p = Point::from_x(x);
            let p = self.mult(&p, &BigUint::from(4u32)); // remove even factor from order
            if p.is_zero() { continue }

            let mut order_divisor = BigUint::one();

            let supersingular = self.is_supersingular_inner(0, PRIMES.len(), &p, &mut order_divisor, &BigUint::from(4u32));
            if supersingular.is_some() {
                return supersingular.unwrap();
            }
            
            // failed to confirm or refute supersingularity
            // try again with another random point
        }
    }
}

// Test MontgomeryCurve against computations verified in SageMath
#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use super::*;

    #[test]
    fn add3() {
        let m = BigUint::from_str("4440482909658554975370291734803070876683357307701517985994904950221288582971213972383728288417044196240794423766370942517241078293092257995419579746810509").unwrap();
        let mut e = MontgomeryCurve::new(BigUint::zero());
        e.a.x = (&e.a.x).mulm(&m, &P);
        e.a.z = (&e.a.z).mulm(&m, &P);

        let p = Point::from_x(BigUint::from_str("2312358140734276205846945336837851905999822779435115388061544163252550153570538668821946365853850479473373914064657421785837590798620029506408147832689588").unwrap());
        let q = Point::from_x(BigUint::from_str("1966862861234634901995729567996927158856250002428122626145328998447532127805166091065319349076607345664792832280295348278073162978962961909004455816564443").unwrap());
        let pq_sum = Point::from_x(BigUint::from_str("4496438020496766616282994673375498894252483353026319191125758947004822107070182520645733369602767683759440544875273015416725059634135174914175766939012166").unwrap());
        let pq_diff = Point::from_x(BigUint::from_str("4731032281376241832924728066192436626155172008235730893583120849924691131260900844876429959545736397402621677725711256678314518853142702763386235211151347").unwrap());

        assert_eq!(e.add3(&p, &q, &pq_sum).normalize(), pq_diff);
        assert_eq!(e.add3(&p, &q, &pq_diff).normalize(), pq_sum);

        let a = BigUint::from_str("3761345407298064040496734252078593163266383672948683311982036640586525413088935630843747460636021803468026326856059472255583832635669585412131689148241997").unwrap();
        let m = BigUint::from_str("4335018931555600401467435874834192213671098660508933828743474302758081805946241114124646029637940922473072839408602410334750139821682401047466066441520451").unwrap();
        let mut e = MontgomeryCurve::new(a);
        e.a.x = (&e.a.x).mulm(&m, &P);
        e.a.z = (&e.a.z).mulm(&m, &P);

        let p = Point::from_x(BigUint::from_str("4495586086716003855534597567300884452536691491748772275003589245631180523337279652783627888301578150162888549581407021085079723193331174354529450141813801").unwrap());
        let q = Point::from_x(BigUint::from_str("2968156693709759315089228567246446056329032951687538433806112532252530237695395332247067250743253776088238540204013357786035386112613018900138010737744091").unwrap());
        let pq_sum = Point::from_x(BigUint::from_str("4836568961178951459879098012443938853977534455317225177432516118051414261523118981615417779011677914332941757424552957519998224019954574117876006773867965").unwrap());
        let pq_diff = Point::from_x(BigUint::from_str("1220043028688943352752850925142751914412058201476466295874059652032179108840941156661901057029693822874756623887252828932021847804313508427021299452110722").unwrap());

        assert_eq!(e.add3(&p, &q, &pq_sum).normalize(), pq_diff);
        assert_eq!(e.add3(&p, &q, &pq_diff).normalize(), pq_sum);
    }

    #[test]
    fn double() {
        let m = BigUint::from_str("521432353287423532365683116793916227093898189955855699219856484683961168662472999598878812772343468478351575748149742789958602618085078626664100234470793").unwrap();
        let mut e = MontgomeryCurve::new(BigUint::zero());
        e.a.x = (&e.a.x).mulm(&m, &P);
        e.a.z = (&e.a.z).mulm(&m, &P);

        let p = Point::from_x(BigUint::from_str("1937602058922554376729573884420059713904515504221308236688786931373958308013090692618208128532656870205905241607945790144785759482324841103257085151564185").unwrap());
        let q = Point::from_x(BigUint::from_str("4216250901005565432547001457763459901579976118362895665057193106411438683120995842641766367782645354045134370935541922935163580753484638504884617740357371").unwrap());
        let two_p = Point::from_x(BigUint::from_str("461982659298801001718203982104366362282888453764857564810571685099923594149330642707265790787585888193988279975777958118774278575028105897433152294744769").unwrap());
        let two_q = Point::from_x(BigUint::from_str("145744890059114125125560197215401329492989650142916373755678630641880230493012904236143199250318243547726234578644878390189910748832599495244142194587720").unwrap());

        assert_eq!(e.double(&p).normalize(), two_p);
        assert_eq!(e.double(&q).normalize(), two_q);

        let a = BigUint::from_str("1709179145736732726190467057772287838888749898111981175668039410052353567995718751503892000081725932287177755925719395202935581874824412236400752119437261").unwrap();
        let m = BigUint::from_str("3232089984278901293276106137927803184522498662311616807660230412533114641533833454660164335995853215609150979443920597833560932719735561181555131755781527").unwrap();
        let mut e = MontgomeryCurve::new(a);
        e.a.x = (&e.a.x).mulm(&m, &P);
        e.a.z = (&e.a.z).mulm(&m, &P);

        let p = Point::from_x(BigUint::from_str("1809210106972348290992138466942598742939032607295145037532561831011318955704255610491397207345121502684825264759784036949158705639565797064768414045286843").unwrap());
        let q = Point::from_x(BigUint::from_str("3316622817677519818985121645043921716174787690029019293216929147225736194796978692153860873593587040902929679904313311019990000675876868194561411288552914").unwrap());
        let two_p = Point::from_x(BigUint::from_str("2157251028612023277034255784142877994688790938730421069519630989106955855484812206272784343211746696219577394029905473040697595489912545789451826952997447").unwrap());
        let two_q = Point::from_x(BigUint::from_str("746635998476565278471879668573951751397281056925627955991523080051482777875006641871676157042375843280167956890761965317435974624098761025237579999741752").unwrap());

        assert_eq!(e.double(&p).normalize(), two_p);
        assert_eq!(e.double(&q).normalize(), two_q);
    }

    #[test]
    fn mult() {
        let a = BigUint::from_str("806328403495992267109149450868966636991037459819762153296867249317594469447447955344137556739687998616446227640966045492819498968091824989742441223431822").unwrap();
        let m = BigUint::from_str("2343990914676338216493302482412120769688102928092223072751560192912394080870946318692017047586271463306051099687255773434603381092195421003414637549482316").unwrap();
        let mut e = MontgomeryCurve::new(a);
        e.a.x = (&e.a.x).mulm(&m, &P);
        e.a.z = (&e.a.z).mulm(&m, &P);

        let p = Point::from_x(BigUint::from_str("1538479638774176874430290160819066779956824411302506566340952445920164206410596012176206212551212863745949091286457738164909520235809167593642776256451348").unwrap());
        let q = Point::from_x(BigUint::from_str("5213933254543299671154526010994617680832632379959926545922182614816102995395142995205368221299305516955759965074840683014659775047482248220893460782049355").unwrap());

        let k = BigUint::zero();
        assert_eq!(e.mult(&p, &k).normalize(), Point::zero());
        assert_eq!(e.mult(&q, &k).normalize(), Point::zero());

        let k = BigUint::one();
        assert_eq!(e.mult(&p, &k).normalize(), p);
        assert_eq!(e.mult(&q, &k).normalize(), q);

        let k = BigUint::from_str("1413548981623491970421649012354808654204587476429295693715967402719052031331206857389616701585512860536287680679987038959402123819391162288336558202497316").unwrap();
        let kp = Point::from_x(BigUint::from_str("2868756308881603928495903678299702212093180571231991285456465524452348559594888946419265166874694035131722161670580533023868466144866644686925596344818041").unwrap());
        let kq = Point::from_x(BigUint::from_str("3405406050475848249398605003979403540220313857436751644846308659849618781571445978170966794736884105229590982900649475807077704413192004717224334412005047").unwrap());
        assert_eq!(e.mult(&p, &k).normalize(), kp);
        assert_eq!(e.mult(&q, &k).normalize(), kq);

        let k = BigUint::from_str("696230786747357566586975788811703550671222779362105476326241069287539761285997871738936119408694768119547302957429166517494205025415697128313440505106924").unwrap();
        let kp = Point::from_x(BigUint::from_str("1065322933146839210816150372813911496132449702212035899660235387485984997901483608944975099026409201713384588143995018893871599603978334571300293868209716").unwrap());
        let kq = Point::from_x(BigUint::from_str("3304617941839255200664447090553691513042451086372604339615061582936491261220222966922855464818834627752367446039290843143082764304175855178935927944912457").unwrap());
        assert_eq!(e.mult(&p, &k).normalize(), kp);
        assert_eq!(e.mult(&q, &k).normalize(), kq);
    }

    #[test]
    fn is_nonsingular() {
        let m = BigUint::from_str("1482991386244245872268337542093393670124232886460077199859913904910083103686769133901993388263105995354924737291146218463428480344871914467085963068133434").unwrap();
        let mut e = MontgomeryCurve::new(BigUint::zero());
        e.a.x = (&e.a.x).mulm(&m, &P);
        e.a.z = (&e.a.z).mulm(&m, &P);
        assert!(e.is_nonsingular());

        let m = BigUint::from_str("3526136202670602223015248024225370523557324558553832722480013465512013502292838617573292446014142512422071366478152209047727336057076022111122658623318445").unwrap();
        let mut e = MontgomeryCurve::new(BigUint::from(2u32));
        e.a.x = (&e.a.x).mulm(&m, &P);
        e.a.z = (&e.a.z).mulm(&m, &P);
        assert!(!e.is_nonsingular());

        let minus_two: BigUint = (&*P) - BigUint::from(2u32);
        let m = BigUint::from_str("4886679815737153389578581995867577711656749121478120131465408385131447251027037051061802473946842453641214615041362243087372292624385854409020419507039262").unwrap();
        let mut e = MontgomeryCurve::new(minus_two);
        e.a.x = (&e.a.x).mulm(&m, &P);
        e.a.z = (&e.a.z).mulm(&m, &P);
        assert!(!e.is_nonsingular());

        let a = BigUint::from_str("806328403495992267109149450868966636991037459819762153296867249317594469447447955344137556739687998616446227640966045492819498968091824989742441223431822").unwrap();
        let m = BigUint::from_str("3524217602009234146288585319860811799255230600204598256003793710987356501758427197666817831234171838927657430163289481096997047220255617499981248071161109").unwrap();
        let mut e = MontgomeryCurve::new(a);
        e.a.x = (&e.a.x).mulm(&m, &P);
        e.a.z = (&e.a.z).mulm(&m, &P);

        assert!(e.is_nonsingular());
        
        let a = BigUint::from_str("1709179145736732726190467057772287838888749898111981175668039410052353567995718751503892000081725932287177755925719395202935581874824412236400752119437261").unwrap();
        let m = BigUint::from_str("2182878533292957361353381006739858886840562083800283831070487915674444875666144618620877871580459773828684117700969188827256220111033815732180495204548636").unwrap();
        let mut e = MontgomeryCurve::new(a);
        e.a.x = (&e.a.x).mulm(&m, &P);
        e.a.z = (&e.a.z).mulm(&m, &P);

        assert!(e.is_nonsingular());
    }

    #[test]
    fn is_supersingular() {
        // Supersingular curves
        let m = BigUint::from_str("4561122714804072297964248928604431772071226190592676084116934519858418127468031603553901854927880992502690163637324040606314250216034085628238683026966664").unwrap();
        let mut e = MontgomeryCurve::new(BigUint::zero());
        e.a.x = (&e.a.x).mulm(&m, &P);
        e.a.z = (&e.a.z).mulm(&m, &P);
        assert!(e.is_supersingular());

        let a = BigUint::from_str("3761345407298064040496734252078593163266383672948683311982036640586525413088935630843747460636021803468026326856059472255583832635669585412131689148241997").unwrap();
        let m = BigUint::from_str("4507816919602490090449930683544672887445469825246471280381825915796162109231566918390173277642686922482234962741570874118413112461758883353012371837338011").unwrap();
        let mut e = MontgomeryCurve::new(a);
        e.a.x = (&e.a.x).mulm(&m, &P);
        e.a.z = (&e.a.z).mulm(&m, &P);
        assert!(e.is_supersingular());

        let a = BigUint::from_str("806328403495992267109149450868966636991037459819762153296867249317594469447447955344137556739687998616446227640966045492819498968091824989742441223431822").unwrap();
        let m = BigUint::from_str("682589759858867006019945921981683911930844413103153096751081265553732323158453661898923474704061121351403714997566638936235758098864420090141210600541605").unwrap();
        let mut e = MontgomeryCurve::new(a);
        e.a.x = (&e.a.x).mulm(&m, &P);
        e.a.z = (&e.a.z).mulm(&m, &P);
        assert!(e.is_supersingular());
        
        let a = BigUint::from_str("1709179145736732726190467057772287838888749898111981175668039410052353567995718751503892000081725932287177755925719395202935581874824412236400752119437261").unwrap();
        let m = BigUint::from_str("1372299461168028207574437529954930118298843176446652665420970385704081105088164666433577075908289791477766172637812347564025927002291465951964304934275749").unwrap();
        let mut e = MontgomeryCurve::new(a);
        e.a.x = (&e.a.x).mulm(&m, &P);
        e.a.z = (&e.a.z).mulm(&m, &P);
        assert!(e.is_supersingular());

        // Ordinary curves
        let a = BigUint::from_str("1939902019534806018427196602142660034992970762099605557756623858435980704369988592447834924216529998385369289786266651499290995707100293528435489458165344").unwrap();
        let m = BigUint::from_str("4043127932910274485564126158378956023446121843039928199364776041672124944895083402139733460354243792206387730324075259812478984755212444548026218178920658").unwrap();
        let mut e = MontgomeryCurve::new(a);
        e.a.x = (&e.a.x).mulm(&m, &P);
        e.a.z = (&e.a.z).mulm(&m, &P);
        assert!(!e.is_supersingular());

        let a = BigUint::from_str("1866223479810078799335388829853753029730709787935267021757606377075016878572548078343651839809404146780592040359635385767112537414085805658738687484339631").unwrap();
        let m = BigUint::from_str("3640326944036693676475257727036086990099851523941743651295162127638831471205355580550155137235405452465381206319509761479787085502435824357185503728478626").unwrap();
        let mut e = MontgomeryCurve::new(a);
        e.a.x = (&e.a.x).mulm(&m, &P);
        e.a.z = (&e.a.z).mulm(&m, &P);
        assert!(!e.is_supersingular());

        let a = BigUint::from_str("2661785568965525579235733734106882226803961945735665715347090419442451011470046364314754336833920717330106734294787860943873523408164640783182492955709999").unwrap();
        let m = BigUint::from_str("683201296837654516155062671807374643518145895323214941072784881241198594642968581707421521295675835749259753224602137584501907664571532242696362741510003").unwrap();
        let mut e = MontgomeryCurve::new(a);
        e.a.x = (&e.a.x).mulm(&m, &P);
        e.a.z = (&e.a.z).mulm(&m, &P);
        assert!(!e.is_supersingular());
    }
}
