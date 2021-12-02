use std::{
    io::{Read, Write, BufReader, BufWriter, BufRead, Seek},
    fs::{File, metadata},
    cmp::{Ordering, min},
    path::Path,
    time::Instant,
    env,
};
    

// Convenience functions for buffered I/O ---------------------------------------------------------- Convenience functions for buffered I/O
#[derive(PartialEq, Eq)]
enum BufferState {
    NotEmpty,
    Empty,
}

trait BufferedRead {
    fn read_byte(&mut self, input: &mut [u8; 1]);
    fn read_usize(&mut self, input: &mut [u8; 8]);
    fn fill_buffer(&mut self) -> BufferState;
}
impl BufferedRead for BufReader<File> {
    fn read_byte(&mut self, input: &mut [u8; 1]) {
        match self.read(input) {
            Ok(_)  => {},
            Err(e) => { 
                println!("Function read_byte failed."); 
                println!("Error: {}", e);
            },
        };
        if self.buffer().len() <= 0 { 
            self.consume(self.capacity()); 
            match self.fill_buf() {
                Ok(_)  => {},
                Err(e) => {
                    println!("Function read_byte failed.");
                    println!("Error: {}", e);
                },
            }
        }
    }
    fn read_usize(&mut self, input: &mut [u8; 8]) {
        match self.read(input) {
            Ok(_)  => {},
            Err(e) => { 
                println!("Function read_usize failed."); 
                println!("Error: {}", e);
            },
        };
        if self.buffer().len() <= 0 { 
            self.consume(self.capacity()); 
            match self.fill_buf() {
                Ok(_)  => {},
                Err(e) => { 
                    println!("Function read_usize failed."); 
                    println!("Error: {}", e);
                },
            }
        }
    }
    fn fill_buffer(&mut self) -> BufferState {
        self.consume(self.capacity());
        match self.fill_buf() {
            Ok(_)  => {},
            Err(e) => { 
                println!("Function fill_buffer failed."); 
                println!("Error: {}", e);
            },
        }
        if self.buffer().is_empty() { 
            return BufferState::Empty; 
        }
        BufferState::NotEmpty
    }
}
trait BufferedWrite {
    fn write_byte(&mut self, output: u8);
    fn write_usize(&mut self, output: usize);
    fn flush_buffer(&mut self);
}
impl BufferedWrite for BufWriter<File> {
    fn write_byte(&mut self, output: u8) {
        match self.write(&[output]) {
            Ok(_)  => {},
            Err(e) => { 
                println!("Function write_byte failed."); 
                println!("Error: {}", e);
            },
        }
        if self.buffer().len() >= self.capacity() { 
            match self.flush() {
                Ok(_)  => {},
                Err(e) => { 
                    println!("Function write_byte failed."); 
                    println!("Error: {}", e);
                },
            } 
        }
    }
    fn write_usize(&mut self, output: usize) {
        match self.write(&output.to_le_bytes()[..]) {
            Ok(_)  => {},
            Err(e) => { 
                println!("Function write_usize failed."); 
                println!("Error: {}", e);
            },
        }
        if self.buffer().len() >= self.capacity() { 
            match self.flush() {
                Ok(_)  => {},
                Err(e) => { 
                    println!("Function write_usize failed."); 
                    println!("Error: {}", e);
                },
            } 
        }
    }
    fn flush_buffer(&mut self) {
        match self.flush() {
            Ok(_)  => {},
            Err(e) => { 
                println!("Function flush_buffer failed."); 
                println!("Error: {}", e);
            },
        }    
    }
}
fn new_input_file(capacity: usize, file_name: &str) -> BufReader<File> {
    BufReader::with_capacity(capacity, File::open(file_name).unwrap())
}
fn new_output_file(capacity: usize, file_name: &str) -> BufWriter<File> {
    BufWriter::with_capacity(capacity, File::create(file_name).unwrap())
}
// ----------------------------------------------------------------------------------------------------------------------------------------


// Burrows-Wheeler Transform ------------------------------------------------------------------------------------ Burrows-Wheeler Transform
fn block_compare(a: usize, b: usize, block: &[u8]) -> Ordering {
    assert!(a < block.len() && b < block.len());

    let min = min(block[a..].len(), block[b..].len());
    
    // Lexicographical comparison
    let result = block[a..a + min].cmp(
                &block[b..b + min]    );
    
    // Implement wraparound if needed
    if result == Ordering::Equal {
        return [&block[a + min..], &block[0..a]].concat().cmp(
              &[&block[b + min..], &block[0..b]].concat()    );
    }
    result   
}

// Takes in a file, outputs a BWT block. The block consists of the 
// transformed data followed by an 8 byte primary index. 
fn bwt_transform_block(file_in: &mut BufReader<File>) -> Vec<u8> { 
    let mut primary_index: usize = 0; // Starting point for inverse transform
    let mut indexes: Vec<u32> = vec![0; file_in.buffer().len()]; // Indexes into block
    let mut bwt: Vec<u8> = vec![0; file_in.buffer().len()];      // BWT output

    // Create indexes into block
    for i in 0..indexes.len() { indexes[i as usize] = i as u32; }
    
    // Sort indexes
    indexes[..].sort_by(|a, b| block_compare(*a as usize, *b as usize, file_in.buffer()));
    
    // Get primary index and BWT output
    for i in 0..bwt.len() {
        if indexes[i] == 1 {
            primary_index = i;
        }
        if indexes[i] == 0 { 
            bwt[i] = file_in.buffer()[file_in.buffer().len() - 1];
        } else {
            bwt[i] = file_in.buffer()[(indexes[i] as usize) - 1];
        }    
    } 

    // Add primary index to end of block
    for i in (0..8).rev() {
        bwt.push(((primary_index >> i*8) & 0xFF) as u8);
    }
    bwt 
}


// Takes in a BWT block, inverts it, and writes 
// the original data to the output file.
fn inverse_bwt_transform(block: &[u8], file_out: &mut BufWriter<File>) {
    let mut transform_vector: Vec<u32> = vec![0; block.len() - 8];

    // Read primary index
    let mut primary_index: usize = 0;
    let mut shift = 0;
    for i in ((block.len() - 8)..(block.len())).rev() {
        primary_index += (block[i] as usize) << shift;
        shift += 8;
    }
    
    let mut counts = [0u32; 256];
    let mut cumul_counts = [0u32; 256];

    // Get number of occurences for each byte
    for i in 0..block.len() - 8 {
        counts[block[i] as usize] += 1;    
    }

    // Get cumulative counts for each byte
    let mut sum = 0;
    for i in 0..256 {
        cumul_counts[i] = sum;
        sum += counts[i];
        counts[i] = 0;
    }

    // Build transformation vector
    for i in 0..block.len() - 8 {
        let index = block[i] as usize; 
        transform_vector[(counts[index] + cumul_counts[index]) as usize] = i as u32;
        counts[index] += 1;
    }
    
    // Invert transform and output original data
    let mut index = primary_index;
    for _ in 0..block.len() - 8 { 
        file_out.write_byte(block[index]);
        index = transform_vector[index] as usize;
    }
    file_out.flush_buffer();
}
// ----------------------------------------------------------------------------------------------------------------------------------------


// Logistic Functions -------------------------------------------------------------------------------------------------- Logistic Functions
// Returns p = 1/(1 + exp(-d))
// d = (-2047..2047), p = (0..4095)
fn squash(mut d: i32) -> i32 {
    const SQUASH_TABLE: [i32; 33] = [
    1,2,3,6,10,16,27,45,73,120,194,310,488,747,1101,
    1546,2047,2549,2994,3348,3607,3785,3901,3975,4022,
    4050,4068,4079,4085,4089,4092,4093,4094];
    if d > 2047  { return 4095; }
    if d < -2047 { return 0;    }
    let w = d & 127;
    d = (d >> 7) + 16;
    (SQUASH_TABLE[d as usize] * (128 - w) + 
     SQUASH_TABLE[(d + 1) as usize] * w + 64) >> 7
}
// Returns p = ln(d/(1-d)) (Inverse of squash)
// d = (0..4095), p = (-2047..2047)
struct Stretch {
    stretch_table: [i16; 4096],
}
impl Stretch {
    fn new() -> Stretch {
        let mut s = Stretch {
            stretch_table: [0; 4096],
        };
        let mut pi = 0;
        for x in -2047..=2047 {
            let i = squash(x);
            for j in pi..=i {
                s.stretch_table[j as usize] = x as i16;
            }
            pi = i + 1;
        }
        s.stretch_table[4095] = 2047;
        s
    }
    fn stretch(&self, d: i32) -> i32 {
        assert!(d < 4096);
        self.stretch_table[d as usize] as i32
    }
}
// ----------------------------------------------------------------------------------------------------------------------------------------


// Adaptive Probability Map -------------------------------------------------------------------------------------- Adaptive Probability Map
struct Apm {
    stretch:    Stretch,
    bin:        usize,    
    num_cxts:   usize, 
    bin_map:    Vec<u16>, // maps each bin to a squashed value
}
impl Apm {
    fn new(n: usize) -> Apm {
        let mut apm = Apm {  
            stretch:    Stretch::new(), 
            bin:        0, // last pr, context
            num_cxts:   n,
            bin_map:    Vec::with_capacity(n * 33),
        };
        apm.bin_map.resize(n * 33, 0);

        for cxt in 0..apm.num_cxts {
            for bin in 0usize..33 {
                apm.bin_map[(cxt * 33) + bin] = if cxt == 0 {
                    (squash(((bin as i32) - 16) * 128) * 16) as u16
                } else {
                    apm.bin_map[bin]
                }
            }
        }
        apm
    }
    fn p(&mut self, bit: i32, rate: i32, mut pr: i32, cxt: usize) -> i32 {
        assert!(bit == 0 || bit == 1 && pr >= 0 && pr < 4096 && cxt < self.num_cxts);
        self.update(bit, rate);
        
        pr = self.stretch.stretch(pr);   // -2047 to 2047
        let interp_wght = pr & 127;      // Interpolation weight (33 points)
        
        // Each context has a corresponding set of 33 bins, and bin is 
        // a specific bin within the set corresponding to the current context
        self.bin = (((pr + 2048) >> 7) + ((cxt as i32) * 33)) as usize;

        (((self.bin_map[self.bin]     as i32) * (128 - interp_wght) ) + 
        ( (self.bin_map[self.bin + 1] as i32) *        interp_wght) ) >> 11
    }
    fn update(&mut self, bit: i32, rate: i32) {
        assert!(bit == 0 || bit == 1 && rate > 0 && rate < 32);
        
        // Variable g controls direction of update (bit = 1: increase, bit = 0: decrease)
        let g: i32 = (bit << 16) + (bit << rate) - bit - bit;
        self.bin_map[self.bin    ] = ( (self.bin_map[self.bin    ] as i32) + 
                                 ((g - (self.bin_map[self.bin    ] as i32)) >> rate) ) as u16;

        self.bin_map[self.bin + 1] = ( (self.bin_map[self.bin + 1] as i32) + 
                                 ((g - (self.bin_map[self.bin + 1] as i32)) >> rate) ) as u16;
    }
}
// ----------------------------------------------------------------------------------------------------------------------------------------


// State Map -------------------------------------------------------------------------------------------------------------------- State Map
const STATE_TABLE: [[u8; 2]; 256] = [
[  1,  2],[  3,  5],[  4,  6],[  7, 10],[  8, 12],[  9, 13],[ 11, 14], // 0
[ 15, 19],[ 16, 23],[ 17, 24],[ 18, 25],[ 20, 27],[ 21, 28],[ 22, 29], // 7
[ 26, 30],[ 31, 33],[ 32, 35],[ 32, 35],[ 32, 35],[ 32, 35],[ 34, 37], // 14
[ 34, 37],[ 34, 37],[ 34, 37],[ 34, 37],[ 34, 37],[ 36, 39],[ 36, 39], // 21
[ 36, 39],[ 36, 39],[ 38, 40],[ 41, 43],[ 42, 45],[ 42, 45],[ 44, 47], // 28
[ 44, 47],[ 46, 49],[ 46, 49],[ 48, 51],[ 48, 51],[ 50, 52],[ 53, 43], // 35
[ 54, 57],[ 54, 57],[ 56, 59],[ 56, 59],[ 58, 61],[ 58, 61],[ 60, 63], // 42
[ 60, 63],[ 62, 65],[ 62, 65],[ 50, 66],[ 67, 55],[ 68, 57],[ 68, 57], // 49
[ 70, 73],[ 70, 73],[ 72, 75],[ 72, 75],[ 74, 77],[ 74, 77],[ 76, 79], // 56
[ 76, 79],[ 62, 81],[ 62, 81],[ 64, 82],[ 83, 69],[ 84, 71],[ 84, 71], // 63
[ 86, 73],[ 86, 73],[ 44, 59],[ 44, 59],[ 58, 61],[ 58, 61],[ 60, 49], // 70
[ 60, 49],[ 76, 89],[ 76, 89],[ 78, 91],[ 78, 91],[ 80, 92],[ 93, 69], // 77
[ 94, 87],[ 94, 87],[ 96, 45],[ 96, 45],[ 48, 99],[ 48, 99],[ 88,101], // 84
[ 88,101],[ 80,102],[103, 69],[104, 87],[104, 87],[106, 57],[106, 57], // 91
[ 62,109],[ 62,109],[ 88,111],[ 88,111],[ 80,112],[113, 85],[114, 87], // 98
[114, 87],[116, 57],[116, 57],[ 62,119],[ 62,119],[ 88,121],[ 88,121], // 105
[ 90,122],[123, 85],[124, 97],[124, 97],[126, 57],[126, 57],[ 62,129], // 112
[ 62,129],[ 98,131],[ 98,131],[ 90,132],[133, 85],[134, 97],[134, 97], // 119
[136, 57],[136, 57],[ 62,139],[ 62,139],[ 98,141],[ 98,141],[ 90,142], // 126
[143, 95],[144, 97],[144, 97],[ 68, 57],[ 68, 57],[ 62, 81],[ 62, 81], // 133
[ 98,147],[ 98,147],[100,148],[149, 95],[150,107],[150,107],[108,151], // 140
[108,151],[100,152],[153, 95],[154,107],[108,155],[100,156],[157, 95], // 147
[158,107],[108,159],[100,160],[161,105],[162,107],[108,163],[110,164], // 154
[165,105],[166,117],[118,167],[110,168],[169,105],[170,117],[118,171], // 161
[110,172],[173,105],[174,117],[118,175],[110,176],[177,105],[178,117], // 168
[118,179],[110,180],[181,115],[182,117],[118,183],[120,184],[185,115], // 175
[186,127],[128,187],[120,188],[189,115],[190,127],[128,191],[120,192], // 182
[193,115],[194,127],[128,195],[120,196],[197,115],[198,127],[128,199], // 189
[120,200],[201,115],[202,127],[128,203],[120,204],[205,115],[206,127], // 196
[128,207],[120,208],[209,125],[210,127],[128,211],[130,212],[213,125], // 203
[214,137],[138,215],[130,216],[217,125],[218,137],[138,219],[130,220], // 210
[221,125],[222,137],[138,223],[130,224],[225,125],[226,137],[138,227], // 217
[130,228],[229,125],[230,137],[138,231],[130,232],[233,125],[234,137], // 224
[138,235],[130,236],[237,125],[238,137],[138,239],[130,240],[241,125], // 231
[242,137],[138,243],[130,244],[245,135],[246,137],[138,247],[140,248], // 238
[249,135],[250, 69],[ 80,251],[140,252],[249,135],[250, 69],[ 80,251], // 245
[140,252],[  0,  0],[  0,  0],[  0,  0]];  // 252

const LIMIT: usize = 127; // Controls rate of adaptation (higher = slower) (0..512)


struct StateMap {
    cxt:           usize,         // Context of last prediction
    cxt_map:       Vec<u32>,      // Maps a context to a prediction and a count
    recipr_table:  [i32; 512],    // Controls the size of each adjustment to cxt_map
}
impl StateMap {
    fn new(n: usize) -> StateMap {
        let mut sm = StateMap { 
            cxt:           0,
            cxt_map:       Vec::with_capacity(n),
            recipr_table:  [0; 512],
        };
        sm.cxt_map.resize(n, 0);

        for pr in sm.cxt_map.iter_mut() {
            *pr = 1 << 31;
        }
        for i in 0..512 { 
            sm.recipr_table[i] = (32_768 / (i + i + 5)) as i32; 
        }
        sm
    }
    fn p(&mut self, bit: i32, cx: usize) -> i32 {
        assert!(bit == 0 || bit == 1);
        self.update(bit);                      // Update prediction for previous context
        self.cxt = cx;
        (self.cxt_map[self.cxt] >> 20) as i32  // Output prediction for new context
    }
    fn update(&mut self, bit: i32) {
        let count: usize = (self.cxt_map[self.cxt] & 511) as usize;  // Low 9 bits
        let prediction: i32 = (self.cxt_map[self.cxt] >> 14) as i32; // High 18 bits

        if count < LIMIT { self.cxt_map[self.cxt] += 1; }

        // Updates cxt_map based on the difference between the predicted and actual bit
        #[allow(overflowing_literals)]
        let high_23_bits: i32 = 0xFFFFFE00;
        self.cxt_map[self.cxt] = self.cxt_map[self.cxt].wrapping_add(
        (((bit << 18) - prediction) * self.recipr_table[count] & high_23_bits) as u32); 
    }
}
// ----------------------------------------------------------------------------------------------------------------------------------------


// Predictor -------------------------------------------------------------------------------------------------------------------- Predictor
struct Predictor {
    cxt:      usize,      apm1: Apm,   
    cxt4:     usize,      apm2: Apm,   
    pr:       i32,        apm3: Apm,   
    state:    [u8; 256],  apm4: Apm,
    sm:       StateMap,   apm5: Apm,   
    run_cxt:  usize,      apm6: Apm,
}
impl Predictor {
    fn new() -> Predictor {
        Predictor {
            cxt:      1,                   apm1: Apm::new(256),
            cxt4:     0,                   apm2: Apm::new(256),
            pr:       2048,                apm3: Apm::new(65536),
            state:    [0; 256],            apm4: Apm::new(256),
            sm:       StateMap::new(256),  apm5: Apm::new(8192), 
            run_cxt:  0,                   apm6: Apm::new(16384),  
        }
    }
    fn p(&mut self) -> i32 { 
        assert!(self.pr >= 0 && self.pr < 4096);
        self.pr 
    } 
    fn update(&mut self, bit: i32) {
        assert!(bit == 0 || bit == 1);
        self.state[self.cxt] = STATE_TABLE[self.state[self.cxt] as usize][bit as usize];

        self.cxt += self.cxt + bit as usize;
        if self.cxt >= 256 {
            self.cxt4 = (self.cxt4 << 8) | (self.cxt - 256);  // Shift new byte into cxt4
            self.cxt = 1;

            // Check if new byte equals previous byte
            if (self.cxt4 ^ (self.cxt4 >> 8)) & 255 == 0 {
                self.run_cxt += 1;
                if self.run_cxt >= 255 {
                    self.run_cxt = 0;
                }
            } else {
                self.run_cxt = 0;
            }       
        }
    
        self.pr = self.sm.p(bit, self.state[self.cxt] as usize);

        self.pr = self.apm1.p(bit, 5, self.pr, self.cxt) +
                  self.apm2.p(bit, 9, self.pr, self.cxt) + 1 >> 1;
        
        self.pr = self.apm3.p(bit, 7, self.pr, self.cxt | (self.cxt4 << 8) & 0xFF00);

        self.pr = self.apm4.p(bit, 8, self.pr, self.run_cxt); 
        
        self.pr = self.apm5.p(bit, 7, self.pr, self.cxt | (self.cxt4 & 0x1F00)) * 3 + self.pr + 2 >> 2;

        self.pr = self.apm6.p(bit, 7, self.pr, 
        ((self.cxt as u32) ^ (((self.cxt4 as u32) & 0xFFFFFF).wrapping_mul(123456791)) >> 18) as usize) 
        + self.pr + 1 >> 1;  
    }
}
// ----------------------------------------------------------------------------------------------------------------------------------------


// Encoder ------------------------------------------------------------------------------------------------------------------------ Encoder
#[derive(PartialEq, Eq)]
enum Mode {
    Compress,
    Decompress,
}

struct Encoder {
    high:       u32,
    low:        u32,
    predictor:  Predictor,
    file_in:    BufReader<File>,
    file_out:   BufWriter<File>,
    x:          u32,
    mode:       Mode,   
}
impl Encoder {
    fn new(file_in: BufReader<File>, file_out: BufWriter<File>, mode: Mode) -> Encoder {
        let mut enc = Encoder {
            high: 0xFFFFFFFF, 
            low: 0, 
            x: 0, 
            predictor: Predictor::new(), 
            file_in, file_out, mode,
        };  
        // Write header placeholder bytes 
        if enc.mode == Mode::Compress {
            enc.file_out.write_usize(0);
            enc.file_out.write_usize(0);
            enc.file_out.write_usize(0);
        }
        enc
    }
    // Compress one bit
    fn compress_bit(&mut self, bit: i32) {
        let mut p = self.predictor.p() as u32;
        if p < 2048 { p += 1; } else {}
        let mid: u32 = self.low + ((self.high - self.low) >> 12) * p 
                       + ((self.high - self.low & 0x0FFF) * p >> 12);
        if bit == 1 {
            self.high = mid;
        } else {
            self.low = mid + 1;
        }
        self.predictor.update(bit);

        while ( (self.high ^ self.low) & 0xFF000000) == 0 {
            self.file_out.write_byte((self.high >> 24) as u8);
            self.high = (self.high << 8) + 255;
            self.low <<= 8;  
        }
    }
    // Decompress and return one bit
    fn decompress_bit(&mut self) -> i32 {
        let mut byte = [0u8; 1];
        let mut p = self.predictor.p() as u32;
        if p < 2048 { p += 1; } else {}
        let mid: u32 = self.low + ((self.high - self.low) >> 12) * p 
                       + ((self.high - self.low & 0x0FFF) * p >> 12);
        let mut bit: i32 = 0;
        if self.x <= mid {
            bit = 1;
            self.high = mid;
        } else {
            self.low = mid + 1;
        }
        self.predictor.update(bit);
        
        while ( (self.high ^ self.low) & 0xFF000000) == 0 {
            self.high = (self.high << 8) + 255;
            self.low <<= 8;
            self.file_in.read_byte(&mut byte); 
            self.x = (self.x << 8) + byte[0] as u32; 
        }
        bit
    }
    fn flush(&mut self) {
        while ( (self.high ^ self.low) & 0xFF000000) == 0 {
            self.file_out.write_byte((self.high >> 24) as u8);
            self.high = (self.high << 8) + 255;
            self.low <<= 8; 
        }
        self.file_out.write_byte((self.high >> 24) as u8);
        self.file_out.flush_buffer();
    }

    fn compress_block(&mut self, block: &[u8]) {
        for byte in block.iter() {
            for i in (0..=7).rev() {
                self.compress_bit(((*byte >> i) & 1) as i32);
            }   
        }
    }
    fn decompress_block(&mut self, block_size: usize) -> Vec<u8> {
        let mut block: Vec<u8> = Vec::with_capacity(block_size); 
        while block.len() < block.capacity() {
            let mut byte: i32 = 1;
            while byte < 256 {
                byte += byte + self.decompress_bit();
            }
            byte -= 256;
            block.push(byte as u8); 
        }
        block
    }

    // Write 24 byte block data header
    fn write_block_data(&mut self, final_block_size: usize, block_size: usize, num_blocks: usize) {
        self.file_out.get_ref().rewind().unwrap();
        self.file_out.write_usize(final_block_size);
        self.file_out.write_usize(block_size);
        self.file_out.write_usize(num_blocks);    
    }
    // Read 24 byte block data header
    fn read_block_data(&mut self) -> (usize, usize, usize) {
        let mut final_block_size = [0u8; 8];
        let mut block_size = [0u8; 8];
        let mut num_blocks = [0u8; 8];
        self.file_in.read_usize(&mut final_block_size);
        self.file_in.read_usize(&mut block_size);
        self.file_in.read_usize(&mut num_blocks);
        (usize::from_le_bytes(final_block_size), 
         usize::from_le_bytes(block_size),
         usize::from_le_bytes(num_blocks)) 
    }

    // Initialize x to first 4 bytes of compressed data
    fn init_x(&mut self) {
        assert!(self.mode == Mode::Decompress);
        for _ in 0..4 {
            let mut byte = [0u8; 1];
            self.file_in.read_byte(&mut byte);
            self.x = (self.x << 8) + byte[0] as u32;
        }
    }
}
// ----------------------------------------------------------------------------------------------------------------------------------------


fn main() {
    let start_time = Instant::now();
    let args: Vec<String> = env::args().collect();

    let e_file_in  = new_input_file(4096, &args[args.len()-2]);
    let e_file_out = new_output_file(4096, &args[args.len()-1]);

    match (&args[1]).as_str() {
        "c" => {  
            let block_size = args[args.len()-3].parse::<usize>().unwrap();
            let mut final_block_size = 0;
            let mut num_blocks = 0;

            let mut file_in = new_input_file(block_size, &args[args.len()-2]);
            let mut enc = Encoder::new(e_file_in, e_file_out, Mode::Compress);

            // Transform and compress -----------------------------------------
            loop {
                if file_in.fill_buffer() == BufferState::Empty { break; }

                let block = bwt_transform_block(&mut file_in);
                final_block_size = block.len();
            
                enc.compress_block(&block);
                num_blocks += 1;
            } 
            // ----------------------------------------------------------------
            enc.flush();
            enc.write_block_data(final_block_size, block_size + 8, num_blocks);
        }
        "d" => {
            let mut file_out = new_output_file(4096, &args[args.len()-1]);
            let mut dec = Encoder::new(e_file_in, e_file_out, Mode::Decompress);

            let (final_block_size, block_size, num_blocks) = dec.read_block_data();
            
            // Run init_x after reading header
            dec.init_x();

            // Decompress and invert transform --------------------------------
            for _ in 0..(num_blocks-1) {
                let block = dec.decompress_block(block_size);
                inverse_bwt_transform(&block, &mut file_out);
            }
            let block = dec.decompress_block(final_block_size);
            inverse_bwt_transform(&block, &mut file_out);
            // ----------------------------------------------------------------

            file_out.flush_buffer();
        }
        _ => { 
            println!("To Compress: c input output");
            println!("To Decompress: d input output");
        }
    } 
    let file_in_size  = metadata(Path::new(&args[args.len()-2])).unwrap().len();
    let file_out_size = metadata(Path::new(&args[args.len()-1])).unwrap().len();
    println!("Finished Compressing.");   
    println!("{} bytes -> {} bytes in {:.2?}", 
    file_in_size, file_out_size, start_time.elapsed()); 
}



