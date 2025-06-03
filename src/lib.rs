mod b_field_element;
use b_field_element::BFieldElement;

mod digest;
use digest::Digest;

mod mds;

mod sponge;
use itertools::Itertools;
pub use sponge::{Domain, Sponge};

use num_traits::{ConstOne, ConstZero};

pub const STATE_SIZE: usize = 16;
pub const NUM_SPLIT_AND_LOOKUP: usize = 4;
pub const RATE: usize = 10;
pub const NUM_ROUNDS: usize = 7;

/// The lookup table with a high algebraic degree used in the TIP-5 permutation. To verify its
/// correctness, see the test “lookup_table_is_correct.”
const LOOKUP_TABLE: [u8; 256] = [
    0, 7, 26, 63, 124, 215, 85, 254, 214, 228, 45, 185, 140, 173, 33, 240, 29, 177, 176, 32, 8,
    110, 87, 202, 204, 99, 150, 106, 230, 14, 235, 128, 213, 239, 212, 138, 23, 130, 208, 6, 44,
    71, 93, 116, 146, 189, 251, 81, 199, 97, 38, 28, 73, 179, 95, 84, 152, 48, 35, 119, 49, 88,
    242, 3, 148, 169, 72, 120, 62, 161, 166, 83, 175, 191, 137, 19, 100, 129, 112, 55, 221, 102,
    218, 61, 151, 237, 68, 164, 17, 147, 46, 234, 203, 216, 22, 141, 65, 57, 123, 12, 244, 54, 219,
    231, 96, 77, 180, 154, 5, 253, 133, 165, 98, 195, 205, 134, 245, 30, 9, 188, 59, 142, 186, 197,
    181, 144, 92, 31, 224, 163, 111, 74, 58, 69, 113, 196, 67, 246, 225, 10, 121, 50, 60, 157, 90,
    122, 2, 250, 101, 75, 178, 159, 24, 36, 201, 11, 243, 132, 198, 190, 114, 233, 39, 52, 21, 209,
    108, 238, 91, 187, 18, 104, 194, 37, 153, 34, 200, 143, 126, 155, 236, 118, 64, 80, 172, 89,
    94, 193, 135, 183, 86, 107, 252, 13, 167, 206, 136, 220, 207, 103, 171, 160, 76, 182, 227, 217,
    158, 56, 174, 4, 66, 109, 139, 162, 184, 211, 249, 47, 125, 232, 117, 43, 16, 42, 127, 20, 241,
    25, 149, 105, 156, 51, 53, 168, 145, 247, 223, 79, 78, 226, 15, 222, 82, 115, 70, 210, 27, 41,
    1, 170, 40, 131, 192, 229, 248, 255,
];

const ROUND_CONSTANTS: [BFieldElement; NUM_ROUNDS * STATE_SIZE] = [
    // 1st round constants
    BFieldElement::new(1332676891236936200),
    BFieldElement::new(16607633045354064669),
    BFieldElement::new(12746538998793080786),
    BFieldElement::new(15240351333789289931),
    BFieldElement::new(10333439796058208418),
    BFieldElement::new(986873372968378050),
    BFieldElement::new(153505017314310505),
    BFieldElement::new(703086547770691416),
    BFieldElement::new(8522628845961587962),
    BFieldElement::new(1727254290898686320),
    BFieldElement::new(199492491401196126),
    BFieldElement::new(2969174933639985366),
    BFieldElement::new(1607536590362293391),
    BFieldElement::new(16971515075282501568),
    BFieldElement::new(15401316942841283351),
    BFieldElement::new(14178982151025681389),
    // 2nd round constants
    BFieldElement::new(2916963588744282587),
    BFieldElement::new(5474267501391258599),
    BFieldElement::new(5350367839445462659),
    BFieldElement::new(7436373192934779388),
    BFieldElement::new(12563531800071493891),
    BFieldElement::new(12265318129758141428),
    BFieldElement::new(6524649031155262053),
    BFieldElement::new(1388069597090660214),
    BFieldElement::new(3049665785814990091),
    BFieldElement::new(5225141380721656276),
    BFieldElement::new(10399487208361035835),
    BFieldElement::new(6576713996114457203),
    BFieldElement::new(12913805829885867278),
    BFieldElement::new(10299910245954679423),
    BFieldElement::new(12980779960345402499),
    BFieldElement::new(593670858850716490),
    // 3rd round constants
    BFieldElement::new(12184128243723146967),
    BFieldElement::new(1315341360419235257),
    BFieldElement::new(9107195871057030023),
    BFieldElement::new(4354141752578294067),
    BFieldElement::new(8824457881527486794),
    BFieldElement::new(14811586928506712910),
    BFieldElement::new(7768837314956434138),
    BFieldElement::new(2807636171572954860),
    BFieldElement::new(9487703495117094125),
    BFieldElement::new(13452575580428891895),
    BFieldElement::new(14689488045617615844),
    BFieldElement::new(16144091782672017853),
    BFieldElement::new(15471922440568867245),
    BFieldElement::new(17295382518415944107),
    BFieldElement::new(15054306047726632486),
    BFieldElement::new(5708955503115886019),
    // 4th round constants
    BFieldElement::new(9596017237020520842),
    BFieldElement::new(16520851172964236909),
    BFieldElement::new(8513472793890943175),
    BFieldElement::new(8503326067026609602),
    BFieldElement::new(9402483918549940854),
    BFieldElement::new(8614816312698982446),
    BFieldElement::new(7744830563717871780),
    BFieldElement::new(14419404818700162041),
    BFieldElement::new(8090742384565069824),
    BFieldElement::new(15547662568163517559),
    BFieldElement::new(17314710073626307254),
    BFieldElement::new(10008393716631058961),
    BFieldElement::new(14480243402290327574),
    BFieldElement::new(13569194973291808551),
    BFieldElement::new(10573516815088946209),
    BFieldElement::new(15120483436559336219),
    // 5th round constants
    BFieldElement::new(3515151310595301563),
    BFieldElement::new(1095382462248757907),
    BFieldElement::new(5323307938514209350),
    BFieldElement::new(14204542692543834582),
    BFieldElement::new(12448773944668684656),
    BFieldElement::new(13967843398310696452),
    BFieldElement::new(14838288394107326806),
    BFieldElement::new(13718313940616442191),
    BFieldElement::new(15032565440414177483),
    BFieldElement::new(13769903572116157488),
    BFieldElement::new(17074377440395071208),
    BFieldElement::new(16931086385239297738),
    BFieldElement::new(8723550055169003617),
    BFieldElement::new(590842605971518043),
    BFieldElement::new(16642348030861036090),
    BFieldElement::new(10708719298241282592),
    // 6th round constants
    BFieldElement::new(12766914315707517909),
    BFieldElement::new(11780889552403245587),
    BFieldElement::new(113183285481780712),
    BFieldElement::new(9019899125655375514),
    BFieldElement::new(3300264967390964820),
    BFieldElement::new(12802381622653377935),
    BFieldElement::new(891063765000023873),
    BFieldElement::new(15939045541699412539),
    BFieldElement::new(3240223189948727743),
    BFieldElement::new(4087221142360949772),
    BFieldElement::new(10980466041788253952),
    BFieldElement::new(18199914337033135244),
    BFieldElement::new(7168108392363190150),
    BFieldElement::new(16860278046098150740),
    BFieldElement::new(13088202265571714855),
    BFieldElement::new(4712275036097525581),
    // 7th round constants
    BFieldElement::new(16338034078141228133),
    BFieldElement::new(1455012125527134274),
    BFieldElement::new(5024057780895012002),
    BFieldElement::new(9289161311673217186),
    BFieldElement::new(9401110072402537104),
    BFieldElement::new(11919498251456187748),
    BFieldElement::new(4173156070774045271),
    BFieldElement::new(15647643457869530627),
    BFieldElement::new(15642078237964257476),
    BFieldElement::new(1405048341078324037),
    BFieldElement::new(3059193199283698832),
    BFieldElement::new(1605012781983592984),
    BFieldElement::new(7134876918849821827),
    BFieldElement::new(5796994175286958720),
    BFieldElement::new(7251651436095127661),
    BFieldElement::new(4565856221886323991),
];

/// The defining, first column of the (circulant) MDS matrix.
/// Derived from the SHA-256 hash of the ASCII string “Tip5” by dividing the digest into 16-bit
/// chunks.
pub const MDS_MATRIX_FIRST_COLUMN: [i64; STATE_SIZE] = [
    61402, 1108, 28750, 33823, 7454, 43244, 53865, 12034, 56951, 27521, 41351, 40901, 12021, 59689,
    26798, 17845,
];

pub struct Tip5 {
    pub state: [BFieldElement; STATE_SIZE],
}

impl Tip5 {
    #[inline]
    pub const fn new(domain: Domain) -> Self {
        use Domain::*;

        let mut state = [BFieldElement::ZERO; STATE_SIZE];

        match domain {
            VariableLength => (),
            FixedLength => {
                let mut i = RATE;
                while i < STATE_SIZE {
                    state[i] = BFieldElement::ONE;
                    i += 1;
                }
            }
        }

        Self { state }
    }

    #[inline]
    pub const fn offset_fermat_cube_map(x: u16) -> u16 {
        let xx = (x + 1) as u64;
        let xxx = xx * xx * xx;
        ((xxx + 256) % 257) as u16
    }

    #[inline]
    fn split_and_lookup(element: &mut BFieldElement) {
        // let value = element.value();
        let mut bytes = element.raw_bytes();

        #[allow(clippy::needless_range_loop)] // faster like so
        for i in 0..8 {
            // bytes[i] = Self::offset_fermat_cube_map(bytes[i].into()) as u8;
            bytes[i] = LOOKUP_TABLE[bytes[i] as usize];
        }

        *element = BFieldElement::from_raw_bytes(&bytes);
    }

    #[inline(always)]
    fn mds_generated(&mut self) {
        let mut lo: [u64; STATE_SIZE] = [0; STATE_SIZE];
        let mut hi: [u64; STATE_SIZE] = [0; STATE_SIZE];
        for i in 0..STATE_SIZE {
            let b = self.state[i].raw_u64();
            hi[i] = b >> 32;
            lo[i] = b & 0xffffffffu64;
        }

        lo = mds::generated_function(&lo);
        hi = mds::generated_function(&hi);

        for r in 0..STATE_SIZE {
            let s = (lo[r] >> 4) as u128 + ((hi[r] as u128) << 28);

            let s_hi = (s >> 64) as u64;
            let s_lo = s as u64;

            let (res, over) = s_lo.overflowing_add(s_hi * 0xffffffffu64);

            self.state[r] =
                BFieldElement::from_raw_u64(if over { res + 0xffffffffu64 } else { res });
        }
    }

    #[inline(always)]
    #[allow(clippy::needless_range_loop)]
    fn sbox_layer(&mut self) {
        for i in 0..NUM_SPLIT_AND_LOOKUP {
            Self::split_and_lookup(&mut self.state[i]);
        }

        for i in NUM_SPLIT_AND_LOOKUP..STATE_SIZE {
            let sq = self.state[i] * self.state[i];
            let qu = sq * sq;
            self.state[i] *= sq * qu;
        }
    }

    #[inline(always)]
    fn round(&mut self, round_index: usize) {
        self.sbox_layer();
        self.mds_generated();
        for i in 0..STATE_SIZE {
            self.state[i] += ROUND_CONSTANTS[round_index * STATE_SIZE + i];
        }
    }

    #[inline(always)]
    fn permutation(&mut self) {
        for i in 0..NUM_ROUNDS {
            self.round(i);
        }
    }

    /// Functionally equivalent to [`permutation`](Self::permutation). Returns the trace of
    /// applying the permutation; that is, the initial state of the sponge as well as its state
    /// after each round.
    pub fn trace(&mut self) -> [[BFieldElement; STATE_SIZE]; 1 + NUM_ROUNDS] {
        let mut trace = [[BFieldElement::ZERO; STATE_SIZE]; 1 + NUM_ROUNDS];

        trace[0] = self.state;
        for i in 0..NUM_ROUNDS {
            self.round(i);
            trace[1 + i] = self.state;
        }

        trace
    }

    /// Hash 10 [`BFieldElement`]s.
    ///
    /// There is no input-padding because the input length is fixed.
    ///
    /// When you want to hash together two [`Digest`]s, use [`Self::hash_pair`]
    /// instead. In some rare cases you do want to hash a fixed-length string
    /// of individual [`BFieldElement`]s, which is why this function is exposed.
    ///
    /// See also: [`Self::hash_pair`], [`Self::hash`], [`Self::hash_varlen`].
    pub fn hash_10(input: &[BFieldElement; 10]) -> [BFieldElement; Digest::LEN] {
        let mut sponge = Self::new(Domain::FixedLength);

        // absorb once
        sponge.state[..10].copy_from_slice(input);

        sponge.permutation();

        // squeeze once
        sponge.state[..Digest::LEN].try_into().unwrap()
    }

    /// Hash two [`Digest`]s together.
    ///
    /// This function is syntax sugar for calling [`Self::hash_10`] on the
    /// concatenation of the digests' values.
    ///
    /// See also: [`Self::hash_10`], [`Self::hash`], [`Self::hash_varlen`].
    pub fn hash_pair(left: Digest, right: Digest) -> Digest {
        let mut sponge = Self::new(Domain::FixedLength);
        sponge.state[..Digest::LEN].copy_from_slice(&left.values());
        sponge.state[Digest::LEN..2 * Digest::LEN].copy_from_slice(&right.values());

        sponge.permutation();

        let digest_values = sponge.state[..Digest::LEN].try_into().unwrap();
        Digest::new(digest_values)
    }

    /// Hash a variable-length sequence of [`BFieldElement`].
    ///
    /// This function pads the input as its length is variable.
    ///
    /// Note that [`Self::hash_varlen`] and [`Self::hash_10`] are different
    /// functions, even when the input to the former, after padding, agrees with
    /// the input to the latter. The difference comes from the initial value of
    /// the capacity-part of the state, which in the case of variable-length
    /// hashing is all-ones but in the case of fixed-length hashing is
    /// all-zeroes.
    ///
    /// Prefer [`Self::hash`] whenever an object is being hashed whose type
    /// implements [`BFieldCodec`]. However, such an object is not always
    /// available, which is why this function is exposed.
    ///
    /// See also: [`Self::hash_10`], [`Self::hash_pair`], [`Self::hash`].
    //
    // - Apply the correct padding
    // - [Sponge::pad_and_absorb_all()]
    // - Read the digest from the resulting state.
    pub fn hash_varlen(input: &[BFieldElement]) -> Digest {
        let mut sponge = Self::init();
        sponge.pad_and_absorb_all(input);
        let produce: [BFieldElement; Digest::LEN] =
            (&sponge.state[..Digest::LEN]).try_into().unwrap();

        Digest::new(produce)
    }
}

impl Sponge for Tip5 {
    const RATE: usize = RATE;

    fn init() -> Self {
        Self::new(Domain::VariableLength)
    }

    fn absorb(&mut self, input: [BFieldElement; RATE]) {
        self.state[..RATE]
            .iter_mut()
            .zip_eq(&input)
            .for_each(|(a, &b)| *a = b);

        self.permutation();
    }

    fn squeeze(&mut self) -> [BFieldElement; RATE] {
        let produce: [BFieldElement; RATE] = (&self.state[..RATE]).try_into().unwrap();
        self.permutation();

        produce
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hash_10() {
        let input = [
            2_977_285_544_793_697_764,
            8_573_079_213_791_329_436,
            14_740_515_030_531_427_526,
            14_389_955_978_682_590_192,
            1_689_728_978_827_025_832,
            2_977_285_544_793_697_764,
            8_573_079_213_791_329_436,
            14_740_515_030_531_427_526,
            14_389_955_978_682_590_192,
            1_689_728_978_827_025_832,
        ]
        .map(BFieldElement::new);

        let expected_output: [u64; 5] = [
            4_284_109_133_012_162_799,
            9_948_087_830_738_081_755,
            1_299_341_039_090_705_558,
            10_318_697_670_389_295_510,
            16_411_665_177_385_553_945,
        ];

        let output = Tip5::hash_10(&input).map(|bfe| bfe.raw_u64());
        assert_eq!(
            output, expected_output,
            "output: {output:?},\nexpected: {expected_output:?}"
        );
    }

    #[test]
    fn test_hash_varlen() {
        let input = [1, 0].map(BFieldElement::new);

        let expected_output: [u64; 5] = [
            1_730_770_831_742_798_981,
            2_676_322_185_709_933_211,
            8_329_210_750_824_781_744,
            16_756_092_452_590_401_876,
            3_547_445_316_740_171_466,
        ];

        let output = Tip5::hash_varlen(&input).values().map(|bfe| bfe.value());
        assert_eq!(output, expected_output,);
    }
}
