/*
  ____            _ _   _ ______  _____ 
 |  _ \          | | \ | |  ____|/ ____|
 | |_) | __ _  __| |  \| | |__  | (___  
 |  _ < / _` |/ _` | . ` |  __|  \___ \ 
 | |_) | (_| | (_| | |\  | |____ ____) |
 |____/ \__,_|\__,_|_| \_|______|_____/ 
                                EMULATOR
*/

#include <bits/stdc++.h>
#include <SDL2/SDL.h>
#define NTH_BIT(x, n) (((x) >> (n)) & 1)
typedef uint8_t  u8;  typedef int8_t  s8;
typedef uint16_t u16; typedef int16_t s16;
typedef uint32_t u32; typedef int32_t s32;
typedef uint64_t u64; typedef int64_t s64;
typedef long     cpu_time_t; // CPU clock cycle count
typedef unsigned cpu_addr_t; // 16-bit memory address

static const char* sdl_error(const char* str) {
	const char* sdl_str = SDL_GetError();
	if (sdl_str && *sdl_str) str = sdl_str;
	return str;
}

class Sound_Queue {
private:
	enum { buf_size = 2048 };
	enum { buf_count = 3 };
	typedef short sample_t;
	sample_t* volatile bufs;
	SDL_sem* volatile free_sem;
	int volatile read_buf;
	int write_buf, write_pos;
	bool sound_open;
	sample_t* buf(int index) {
		assert((unsigned)index < buf_count);
		return bufs + (long)index * buf_size;
	};
	void fill_buffer(Uint8* out, int count) {
		if (SDL_SemValue(free_sem) < buf_count - 1) {
			memcpy(out, buf(read_buf), count);
			read_buf = (read_buf + 1) % buf_count;
			SDL_SemPost(free_sem);
		}
		else memset(out, 0, count);
	}
	static void fill_buffer_(void* user_data, Uint8* out, int count) { ((Sound_Queue*) user_data)->fill_buffer(out, count); }
public:
	Sound_Queue() {
		bufs = NULL;
		free_sem = NULL;
		write_buf = 0;
		write_pos = 0;
		read_buf = 0;
		sound_open = false;
	}
	~Sound_Queue() {
		if (sound_open) SDL_PauseAudio(1), SDL_CloseAudio();
		if (free_sem) SDL_DestroySemaphore(free_sem);
		delete [] bufs;
	}	
	const char* init(long sample_rate, int chan_count = 1) { // Initialize with specified sample rate and channel count
		assert(!bufs); // Can only be initialized once
		bufs = new sample_t [(long) buf_size * buf_count];
		if (!bufs) return "Out of memory";
		free_sem = SDL_CreateSemaphore( buf_count - 1 );
		if (!free_sem) return sdl_error("Couldn't create semaphore");
		SDL_AudioSpec as;
		as.freq = sample_rate;
		as.format = AUDIO_S16SYS;
		as.channels = chan_count;
		as.silence = 0;
		as.samples = buf_size;
		as.size = 0;
		as.callback = fill_buffer_;
		as.userdata = this;
		if (SDL_OpenAudio(&as, NULL ) < 0) return sdl_error("Couldn't open SDL audio");
		SDL_PauseAudio(false);
		sound_open = true;
		return NULL;
	}
	int sample_count() const; // Number of samples in buffer waiting to be played
	void write(const sample_t* in, int count) { // Write samples to buffer and block until enough space is available
		while (count) {
			int n = buf_size - write_pos;
			if (n > count) n = count;
			memcpy(buf(write_buf) + write_pos, in, n * sizeof(sample_t));
			in += n, write_pos += n, count -= n;
			if (write_pos >= buf_size) {
				write_pos = 0, write_buf = (write_buf + 1) % buf_count;
				SDL_SemWait(free_sem);
			}
		}
	}
};

struct apu_snapshot_t;
class Nonlinear_Buffer;

class Nes_Apu {
private:
	friend class Nes_Nonlinearizer;
	void enable_nonlinear(double volume);
	Nes_Apu(const Nes_Apu&);
	Nes_Apu& operator = (const Nes_Apu&);
	Nes_Osc* oscs [osc_count];
	Nes_Square square1, square2;
	Nes_Noise noise;
	Nes_Triangle triangle;
	Nes_Dmc dmc;
	cpu_time_t last_time, earliest_irq_, next_irq;
	int frame_period, frame_delay, frame, osc_enables, frame_mode;
	bool irq_flag;
	void (*irq_notifier_)( void* user_data );
	void* irq_data;
	Nes_Square::Synth square_synth;	
	void irq_changed(), state_restored();
	friend struct Nes_Dmc;
public:
	Nes_Apu();
	~Nes_Apu();
	void output(Blip_Buffer*); // Set buffer to generate all sound into, or disable sound if NULL
	inline void dmc_reader(int(*func)(void*, cpu_addr_t), void* user_data = NULL) { // Set memory reader callback used by DMC oscillator to fetch samples
		dmc.rom_reader_data = user_data, dmc.rom_reader = func;
	}
	// Write to register (0x4000-0x4017, except 0x4014 and 0x4016)
	enum { start_addr = 0x4000 };
	enum { end_addr   = 0x4017 };
	void write_register( cpu_time_t, cpu_addr_t, int data );
	// Read from status register at 0x4015
	enum { status_addr = 0x4015 };
	int read_status(cpu_time_t);
	void end_frame(cpu_time_t); // Run all oscillators up to specified time
	void reset( bool pal_timing = false, int initial_dmc_dac = 0 ); // Reset internal frame counter, registers, and all oscillators
	// Save/load snapshot of exact emulation state
	void save_snapshot(apu_snapshot_t* out) const;
	void load_snapshot(apu_snapshot_t const&);
	void buffer_cleared(); // Reset oscillator amplitudes
	void treble_eq(const blip_eq_t&); // Set treble equalization
	// Set sound output of specific oscillator to buffer
	enum { osc_count = 5 };
	void osc_output(int index, Blip_Buffer* buffer);
	void irq_notifier( void (*func)( void* user_data ), void* user_data ) { // Set IRQ time callback
		irq_notifier_ = func, irq_data = user_data;
	}
	// Get time that APU-generated IRQ will occur
	enum { no_irq = LONG_MAX / 2 + 1 };
	enum { irq_waiting = 0 };
	cpu_time_t earliest_irq() const { return earliest_irq_; }
	inline int count_dmc_reads( cpu_time_t t, cpu_time_t* last_read = NULL ) const { return dmc.count_reads( time, last_read ); }
	void run_until( cpu_time_t ); // Run APU until specified time
	inline void osc_output( int osc, Blip_Buffer* buf ) {
		assert(("Nes_Apu::osc_output(): Index out of range", 0 <= osc && osc < osc_count));
		oscs[osc]->output = buf;
	}
};

struct Nes_Osc
{
	unsigned char regs [4];
	bool reg_written [4];
	Blip_Buffer* output;
	int length_counter;// length counter (0 if unused by oscillator)
	int delay;      // delay until next (potential) transition
	int last_amp;   // last amplitude oscillator was outputting
	
	void clock_length( int halt_mask );
	int period() const {
		return (regs [3] & 7) * 0x100 + (regs [2] & 0xff);
	}
	void reset() {
		delay = 0;
		last_amp = 0;
	}
	int update_amp( int amp ) {
		int delta = amp - last_amp;
		last_amp = amp;
		return delta;
	}
};

struct Nes_Envelope : Nes_Osc
{
	int envelope;
	int env_delay;
	
	void clock_envelope();
	int volume() const;
	void reset() {
		envelope = 0;
		env_delay = 0;
		Nes_Osc::reset();
	}
};

// Nes_Square
struct Nes_Square : Nes_Envelope
{
	enum { negate_flag = 0x08 };
	enum { shift_mask = 0x07 };
	enum { phase_range = 8 };
	int phase;
	int sweep_delay;
	
	typedef Blip_Synth<blip_good_quality,15> Synth;
	const Synth* synth; // shared between squares
	
	void clock_sweep( int adjust );
	void run( cpu_time_t, cpu_time_t );
	void reset() {
		sweep_delay = 0;
		Nes_Envelope::reset();
	}
};

// Nes_Triangle
struct Nes_Triangle : Nes_Osc
{
	enum { phase_range = 16 };
	int phase;
	int linear_counter;
	Blip_Synth<blip_good_quality,15> synth;
	
	int calc_amp() const;
	void run( cpu_time_t, cpu_time_t );
	void clock_linear_counter();
	void reset() {
		linear_counter = 0;
		phase = phase_range;
		Nes_Osc::reset();
	}
};

// Nes_Noise
struct Nes_Noise : Nes_Envelope
{
	int noise;
	Blip_Synth<blip_med_quality,15> synth;
	
	void run( cpu_time_t, cpu_time_t );
	void reset() {
		noise = 1 << 14;
		Nes_Envelope::reset();
	}
};

// Nes_Dmc
struct Nes_Dmc : Nes_Osc
{
	int address;    // address of next byte to read
	int period;
	//int length_counter; // bytes remaining to play (already defined in Nes_Osc)
	int buf;
	int bits_remain;
	int bits;
	bool buf_empty;
	bool silence;
	
	enum { loop_flag = 0x40 };
	
	int dac;
	
	cpu_time_t next_irq;
	bool irq_enabled;
	bool irq_flag;
	bool pal_mode;
	bool nonlinear;
	
	int (*rom_reader)( void*, cpu_addr_t ); // needs to be initialized to rom read function
	void* rom_reader_data;
	
	Nes_Apu* apu;
	
	Blip_Synth<blip_med_quality,127> synth;
	
	void start();
	void write_register( int, int );
	void run( cpu_time_t, cpu_time_t );
	void recalc_irq();
	void fill_buffer();
	void reload_sample();
	void reset();
	int count_reads( cpu_time_t, cpu_time_t* ) const;
};

class Blip_Reader;

// Source time unit.
typedef long blip_time_t;

// Type of sample produced. Signed 16-bit format.
typedef int16_t blip_sample_t;

// Make buffer as large as possible (currently about 65000 samples)
const int blip_default_length = 0;

class Blip_Buffer {
public:
	// Construct an empty buffer.
	Blip_Buffer();
	~Blip_Buffer();
	
	// Set output sample rate and buffer length in milliseconds (1/1000 sec),
	// then clear buffer. If length is not specified, make as large as possible.
	// If there is insufficient memory for the buffer, sets the buffer length
	// to 0 and returns error string (or propagates exception if compiler supports it).
	blargg_err_t sample_rate( long samples_per_sec, int msec_length = blip_default_length );
	// to do: rename to set_sample_rate
	
	// Length of buffer, in milliseconds
	int length() const;
	
	// Current output sample rate
	long sample_rate() const;
	
	// Number of source time units per second
	void clock_rate( long );
	long clock_rate() const;
	
	// Set frequency at which high-pass filter attenuation passes -3dB
	void bass_freq( int frequency );
	
	// Remove all available samples and clear buffer to silence. If 'entire_buffer' is
	// false, just clear out any samples waiting rather than the entire buffer.
	void clear( bool entire_buffer = true );
	
	// to do:
	// Notify Blip_Buffer that synthesis has been performed until specified time
	//void run_until( blip_time_t );
	
	// End current time frame of specified duration and make its samples available
	// (along with any still-unread samples) for reading with read_samples(). Begin
	// a new time frame at the end of the current frame. All transitions must have
	// been added before 'time'.
	void end_frame( blip_time_t time );
	
	// Number of samples available for reading with read_samples()
	long samples_avail() const;
	
	// Read at most 'max_samples' out of buffer into 'dest', removing them from from
	// the buffer. Return number of samples actually read and removed. If stereo is
	// true, increment 'dest' one extra time after writing each sample, to allow
	// easy interleving of two channels into a stereo output buffer.
	long read_samples( blip_sample_t* dest, long max_samples, bool stereo = false );
	
	// Remove 'count' samples from those waiting to be read
	void remove_samples( long count );
	
	// Number of samples delay from synthesis to samples read out
	int output_latency() const;
	
	
	// Experimental external buffer mixing support
	
	// Number of raw samples that can be mixed within frame of specified duration
	long count_samples( blip_time_t duration ) const;
	
	// Mix 'count' samples from 'buf' into buffer.
	void mix_samples( const blip_sample_t* buf, long count );
	
	
	// not documented yet
	
	void remove_silence( long count );
	
	typedef unsigned long resampled_time_t;
	
	resampled_time_t resampled_time( blip_time_t t ) const {
		return t * resampled_time_t (factor_) + offset_;
	}
	
	resampled_time_t resampled_duration( int t ) const {
		return t * resampled_time_t (factor_);
	}
	
private:
	// noncopyable
	Blip_Buffer( const Blip_Buffer& );
	Blip_Buffer& operator = ( const Blip_Buffer& );

	// Don't use the following members. They are public only for technical reasons.
	public:
		enum { widest_impulse_ = 24 };
		typedef BOOST::uint16_t buf_t_;
		
		unsigned long factor_;
		resampled_time_t offset_;
		buf_t_* buffer_;
		unsigned buffer_size_;
	private:
		long reader_accum;
		int bass_shift;
		long samples_per_sec;
		long clocks_per_sec;
		int bass_freq_;
		int length_;
		
		enum { accum_fract = 15 }; // less than 16 to give extra sample range
		enum { sample_offset = 0x7F7F }; // repeated byte allows memset to clear buffer
		
		friend class Blip_Reader;
};

// Low-pass equalization parameters (see notes.txt)
class blip_eq_t {
public:
	blip_eq_t( double treble = 0 );
	blip_eq_t( double treble, long cutoff, long sample_rate );
private:
	double treble;
	long cutoff;
	long sample_rate;
	friend class Blip_Impulse_;
};

// not documented yet (see Multi_Buffer.cpp for an example of use)
class Blip_Reader {
	const Blip_Buffer::buf_t_* buf;
	long accum;
	#ifdef __MWERKS__
	void operator = ( struct foobar ); // helps optimizer
	#endif
public:
	// avoid anything which might cause optimizer to put object in memory
	
	int begin( Blip_Buffer& blip_buf ) {
		buf = blip_buf.buffer_;
		accum = blip_buf.reader_accum;
		return blip_buf.bass_shift;
	}
	
	int read() const {
		return accum >> Blip_Buffer::accum_fract;
	}
	
	void next( int bass_shift = 9 ) {
		accum -= accum >> bass_shift;
		accum += ((long) *buf++ - Blip_Buffer::sample_offset) << Blip_Buffer::accum_fract;
	}
	
	void end( Blip_Buffer& blip_buf ) {
		blip_buf.reader_accum = accum;
	}
};



// End of public interface
	
#ifndef BLIP_BUFFER_ACCURACY
	#define BLIP_BUFFER_ACCURACY 16
#endif

const int blip_res_bits_ = 5;

typedef uint32_t blip_pair_t_;

class Blip_Impulse_ {
	typedef uint16_t imp_t;
	
	blip_eq_t eq;
	double  volume_unit_;
	imp_t*  impulses;
	imp_t*  impulse;
	int     width;
	int     fine_bits;
	int     res;
	bool    generate;
	
	void fine_volume_unit();
	void scale_impulse( int unit, imp_t* ) const;
public:
	Blip_Buffer*    buf;
	uint32_t offset;
	
	void init( blip_pair_t_* impulses, int width, int res, int fine_bits = 0 );
	void volume_unit( double );
	void treble_eq( const blip_eq_t& );
};

inline blip_eq_t::blip_eq_t( double t ) :
		treble( t ), cutoff( 0 ), sample_rate( 44100 ) {
}

inline blip_eq_t::blip_eq_t( double t, long c, long sr ) :
		treble( t ), cutoff( c ), sample_rate( sr ) {
}

inline int Blip_Buffer::length() const {
	return length_;
}

inline long Blip_Buffer::samples_avail() const {
	return long (offset_ >> BLIP_BUFFER_ACCURACY);
}

inline long Blip_Buffer::sample_rate() const {
	return samples_per_sec;
}

inline void Blip_Buffer::end_frame( blip_time_t t ) {
	offset_ += t * factor_;
	assert(( "Blip_Buffer::end_frame(): Frame went past end of buffer",
			samples_avail() <= (long) buffer_size_ ));
}

inline void Blip_Buffer::remove_silence( long count ) {
	assert(( "Blip_Buffer::remove_silence(): Tried to remove more samples than available",
			count <= samples_avail() ));
	offset_ -= resampled_time_t (count) << BLIP_BUFFER_ACCURACY;
}

inline int Blip_Buffer::output_latency() const {
	return widest_impulse_ / 2;
}

inline long Blip_Buffer::clock_rate() const {
	return clocks_per_sec;
}

Nes_Apu::Nes_Apu()
{
	dmc.apu = this;
	dmc.rom_reader = NULL;
	square1.synth = &square_synth;
	square2.synth = &square_synth;
	irq_notifier_ = NULL;
	
	oscs [0] = &square1;
	oscs [1] = &square2;
	oscs [2] = &triangle;
	oscs [3] = &noise;
	oscs [4] = &dmc;
	
	output( NULL );
	volume( 1.0 );
	reset( false );
}

Nes_Apu::~Nes_Apu()
{
}

void Nes_Apu::treble_eq( const blip_eq_t& eq )
{
	square_synth.treble_eq( eq );
	triangle.synth.treble_eq( eq );
	noise.synth.treble_eq( eq );
	dmc.synth.treble_eq( eq );
}

void Nes_Apu::buffer_cleared()
{
	square1.last_amp = 0;
	square2.last_amp = 0;
	triangle.last_amp = 0;
	noise.last_amp = 0;
	dmc.last_amp = 0;
}

void Nes_Apu::enable_nonlinear( double v )
{
	dmc.nonlinear = true;
	square_synth.volume( 1.3 * 0.25751258 / 0.742467605 * 0.25 * v );
	
	const double tnd = 0.75 / 202 * 0.48;
	triangle.synth.volume_unit( 3 * tnd );
	noise.synth.volume_unit( 2 * tnd );
	dmc.synth.volume_unit( tnd );
	
	buffer_cleared();
}

void Nes_Apu::volume( double v )
{
	dmc.nonlinear = false;
	square_synth.volume( 0.1128 * v );
	triangle.synth.volume( 0.12765 * v );
	noise.synth.volume( 0.0741 * v );
	dmc.synth.volume( 0.42545 * v );
}

void Nes_Apu::output( Blip_Buffer* buffer )
{
	for ( int i = 0; i < osc_count; i++ )
		osc_output( i, buffer );
}

void Nes_Apu::reset( bool pal_mode, int initial_dmc_dac )
{
	// to do: time pal frame periods exactly
	frame_period = pal_mode ? 8314 : 7458;
	dmc.pal_mode = pal_mode;
	
	square1.reset();
	square2.reset();
	triangle.reset();
	noise.reset();
	dmc.reset();
	
	last_time = 0;
	osc_enables = 0;
	irq_flag = false;
	earliest_irq_ = no_irq;
	frame_delay = 1;
	write_register( 0, 0x4017, 0x00 );
	write_register( 0, 0x4015, 0x00 );
	
	for ( cpu_addr_t addr = start_addr; addr <= 0x4013; addr++ )
		write_register( 0, addr, (addr & 3) ? 0x00 : 0x10 );
	
	dmc.dac = initial_dmc_dac;
	if ( !dmc.nonlinear )
		dmc.last_amp = initial_dmc_dac; // prevent output transition
}

void Nes_Apu::irq_changed()
{
	cpu_time_t new_irq = dmc.next_irq;
	if ( dmc.irq_flag | irq_flag ) {
		new_irq = 0;
	}
	else if ( new_irq > next_irq ) {
		new_irq = next_irq;
	}
	
	if ( new_irq != earliest_irq_ ) {
		earliest_irq_ = new_irq;
		if ( irq_notifier_ )
			irq_notifier_( irq_data );
	}
}

// frames

void Nes_Apu::run_until( cpu_time_t end_time )
{
	require( end_time >= last_time );
	
	if ( end_time == last_time )
		return;
	
	while ( true )
	{
		// earlier of next frame time or end time
		cpu_time_t time = last_time + frame_delay;
		if ( time > end_time )
			time = end_time;
		frame_delay -= time - last_time;
		
		// run oscs to present
		square1.run( last_time, time );
		square2.run( last_time, time );
		triangle.run( last_time, time );
		noise.run( last_time, time );
		dmc.run( last_time, time );
		last_time = time;
		
		if ( time == end_time )
			break; // no more frames to run
		
		// take frame-specific actions
		frame_delay = frame_period;
		switch ( frame++ )
		{
			case 0:
				if ( !(frame_mode & 0xc0) ) {
		 			next_irq = time + frame_period * 4 + 1;
		 			irq_flag = true;
		 		}
		 		// fall through
		 	case 2:
		 		// clock length and sweep on frames 0 and 2
				square1.clock_length( 0x20 );
				square2.clock_length( 0x20 );
				noise.clock_length( 0x20 );
				triangle.clock_length( 0x80 ); // different bit for halt flag on triangle
				
				square1.clock_sweep( -1 );
				square2.clock_sweep( 0 );
		 		break;
		 	
			case 1:
				// frame 1 is slightly shorter
				frame_delay -= 2;
				break;
			
		 	case 3:
		 		frame = 0;
		 		
		 		// frame 3 is almost twice as long in mode 1
		 		if ( frame_mode & 0x80 )
					frame_delay += frame_period - 6;
				break;
		}
		
		// clock envelopes and linear counter every frame
		triangle.clock_linear_counter();
		square1.clock_envelope();
		square2.clock_envelope();
		noise.clock_envelope();
	}
}

void Nes_Apu::end_frame( cpu_time_t end_time )
{
	if ( end_time > last_time )
		run_until( end_time );
	
	// make times relative to new frame
	last_time -= end_time;
	require( last_time >= 0 );
	
	if ( next_irq != no_irq ) {
		next_irq -= end_time;
		assert( next_irq >= 0 );
	}
	if ( dmc.next_irq != no_irq ) {
		dmc.next_irq -= end_time;
		assert( dmc.next_irq >= 0 );
	}
	if ( earliest_irq_ != no_irq ) {
		earliest_irq_ -= end_time;
		if ( earliest_irq_ < 0 )
			earliest_irq_ = 0;
	}
}

// registers

static const unsigned char length_table [0x20] = {
	0x0A, 0xFE, 0x14, 0x02, 0x28, 0x04, 0x50, 0x06,
	0xA0, 0x08, 0x3C, 0x0A, 0x0E, 0x0C, 0x1A, 0x0E, 
	0x0C, 0x10, 0x18, 0x12, 0x30, 0x14, 0x60, 0x16,
	0xC0, 0x18, 0x48, 0x1A, 0x10, 0x1C, 0x20, 0x1E
};

void Nes_Apu::write_register( cpu_time_t time, cpu_addr_t addr, int data )
{
	require( addr > 0x20 ); // addr must be actual address (i.e. 0x40xx)
	require( (unsigned) data <= 0xff );
	
	// Ignore addresses outside range
	if ( addr < start_addr || end_addr < addr )
		return;
	
	run_until( time );
	
	if ( addr < 0x4014 )
	{
		// Write to channel
		int osc_index = (addr - start_addr) >> 2;
		Nes_Osc* osc = oscs [osc_index];
		
		int reg = addr & 3;
		osc->regs [reg] = data;
		osc->reg_written [reg] = true;
		
		if ( osc_index == 4 )
		{
			// handle DMC specially
			dmc.write_register( reg, data );
		}
		else if ( reg == 3 )
		{
			// load length counter
			if ( (osc_enables >> osc_index) & 1 )
				osc->length_counter = length_table [(data >> 3) & 0x1f];
			
			// reset square phase
			if ( osc_index < 2 )
				((Nes_Square*) osc)->phase = Nes_Square::phase_range - 1;
		}
	}
	else if ( addr == 0x4015 )
	{
		// Channel enables
		for ( int i = osc_count; i--; )
			if ( !((data >> i) & 1) )
				oscs [i]->length_counter = 0;
		
		bool recalc_irq = dmc.irq_flag;
		dmc.irq_flag = false;
		
		int old_enables = osc_enables;
		osc_enables = data;
		if ( !(data & 0x10) ) {
			dmc.next_irq = no_irq;
			recalc_irq = true;
		}
		else if ( !(old_enables & 0x10) ) {
			dmc.start(); // dmc just enabled
		}
		
		if ( recalc_irq )
			irq_changed();
	}
	else if ( addr == 0x4017 )
	{
		// Frame mode
		frame_mode = data;
		
		bool irq_enabled = !(data & 0x40);
		irq_flag &= irq_enabled;
		next_irq = no_irq;
		
		// mode 1
		frame_delay = (frame_delay & 1);
		frame = 0;
		
		if ( !(data & 0x80) )
		{
			// mode 0
			frame = 1;
			frame_delay += frame_period;
			if ( irq_enabled )
				next_irq = time + frame_delay + frame_period * 3;
		}
		
		irq_changed();
	}
}

int Nes_Apu::read_status( cpu_time_t time )
{
	run_until( time - 1 );
	
	int result = (dmc.irq_flag << 7) | (irq_flag << 6);
	
	for ( int i = 0; i < osc_count; i++ )
		if ( oscs [i]->length_counter )
			result |= 1 << i;
	
	run_until( time );
	
	if ( irq_flag ) {
		irq_flag = false;
		irq_changed();
	}
	
	return result;
}

void Nes_Osc::clock_length( int halt_mask )
{
	if ( length_counter && !(regs [0] & halt_mask) )
		length_counter--;
}

void Nes_Envelope::clock_envelope()
{
	int period = regs [0] & 15;
	if ( reg_written [3] ) {
		reg_written [3] = false;
		env_delay = period;
		envelope = 15;
	}
	else if ( --env_delay < 0 ) {
		env_delay = period;
		if ( envelope | (regs [0] & 0x20) )
			envelope = (envelope - 1) & 15;
	}
}

int Nes_Envelope::volume() const
{
	return length_counter == 0 ? 0 : (regs [0] & 0x10) ? (regs [0] & 15) : envelope;
}

// Nes_Square

void Nes_Square::clock_sweep( int negative_adjust )
{
	int sweep = regs [1];
	
	if ( --sweep_delay < 0 )
	{
		reg_written [1] = true;
		
		int period = this->period();
		int shift = sweep & shift_mask;
		if ( shift && (sweep & 0x80) && period >= 8 )
		{
			int offset = period >> shift;
			
			if ( sweep & negate_flag )
				offset = negative_adjust - offset;
			
			if ( period + offset < 0x800 )
			{
				period += offset;
				// rewrite period
				regs [2] = period & 0xff;
				regs [3] = (regs [3] & ~7) | ((period >> 8) & 7);
			}
		}
	}
	
	if ( reg_written [1] ) {
		reg_written [1] = false;
		sweep_delay = (sweep >> 4) & 7;
	}
}

void Nes_Square::run( cpu_time_t time, cpu_time_t end_time )
{
	if ( !output )
		return;
	
	const int volume = this->volume();
	const int period = this->period();
	int offset = period >> (regs [1] & shift_mask);
	if ( regs [1] & negate_flag )
		offset = 0;
	
	const int timer_period = (period + 1) * 2;
	if ( volume == 0 || period < 8 || (period + offset) >= 0x800 )
	{
		if ( last_amp ) {
			synth->offset( time, -last_amp, output );
			last_amp = 0;
		}
		
		time += delay;
		if ( time < end_time )
		{
			// maintain proper phase
			int count = (end_time - time + timer_period - 1) / timer_period;
			phase = (phase + count) & (phase_range - 1);
			time += (long) count * timer_period;
		}
	}
	else
	{
		// handle duty select
		int duty_select = (regs [0] >> 6) & 3;
		int duty = 1 << duty_select; // 1, 2, 4, 2
		int amp = 0;
		if ( duty_select == 3 ) {
			duty = 2; // negated 25%
			amp = volume;
		}
		if ( phase < duty )
			amp ^= volume;
		
		int delta = update_amp( amp );
		if ( delta )
			synth->offset( time, delta, output );
		
		time += delay;
		if ( time < end_time )
		{
			Blip_Buffer* const output = this->output;
			const Synth* synth = this->synth;
			int delta = amp * 2 - volume;
			int phase = this->phase;
			
			do {
				phase = (phase + 1) & (phase_range - 1);
				if ( phase == 0 || phase == duty ) {
					delta = -delta;
					synth->offset_inline( time, delta, output );
				}
				time += timer_period;
			}
			while ( time < end_time );
			
			last_amp = (delta + volume) >> 1;
			this->phase = phase;
		}
	}
	
	delay = time - end_time;
}

// Nes_Triangle

void Nes_Triangle::clock_linear_counter()
{
	if ( reg_written [3] )
		linear_counter = regs [0] & 0x7f;
	else if ( linear_counter )
		linear_counter--;
	
	if ( !(regs [0] & 0x80) )
		reg_written [3] = false;
}

inline int Nes_Triangle::calc_amp() const
{
	int amp = phase_range - phase;
	if ( amp < 0 )
		amp = phase - (phase_range + 1);
	return amp;
}

void Nes_Triangle::run( cpu_time_t time, cpu_time_t end_time )
{
	if ( !output )
		return;
	
	// to do: track phase when period < 3
	// to do: Output 7.5 on dac when period < 2? More accurate, but results in more clicks.
	
	int delta = update_amp( calc_amp() );
	if ( delta )
		synth.offset( time, delta, output );
	
	time += delay;
	const int timer_period = period() + 1;
	if ( length_counter == 0 || linear_counter == 0 || timer_period < 3 )
	{
		time = end_time;
	}
	else if ( time < end_time )
	{
		Blip_Buffer* const output = this->output;
		
		int phase = this->phase;
		int volume = 1;
		if ( phase > phase_range ) {
			phase -= phase_range;
			volume = -volume;
		}
		
		do {
			if ( --phase == 0 ) {
				phase = phase_range;
				volume = -volume;
			}
			else {
				synth.offset_inline( time, volume, output );
			}
			
			time += timer_period;
		}
		while ( time < end_time );
		
		if ( volume < 0 )
			phase += phase_range;
		this->phase = phase;
		last_amp = calc_amp();
 	}
	delay = time - end_time;
}

// Nes_Dmc

void Nes_Dmc::reset()
{
	address = 0;
	dac = 0;
	buf = 0;
	bits_remain = 1;
	bits = 0;
	buf_empty = true;
	silence = true;
	next_irq = Nes_Apu::no_irq;
	irq_flag = false;
	irq_enabled = false;
	
	Nes_Osc::reset();
	period = 0x036;
}

void Nes_Dmc::recalc_irq()
{
	cpu_time_t irq = Nes_Apu::no_irq;
	if ( irq_enabled && length_counter )
		irq = apu->last_time + delay +
				((length_counter - 1) * 8 + bits_remain - 1) * cpu_time_t (period) + 1;
	if ( irq != next_irq ) {
		next_irq = irq;
		apu->irq_changed();
	}
}

int Nes_Dmc::count_reads( cpu_time_t time, cpu_time_t* last_read ) const
{
	if ( last_read )
		*last_read = time;
	
	if ( length_counter == 0 )
		return 0; // not reading
	
	long first_read = apu->last_time + delay + long (bits_remain - 1) * period;
	long avail = time - first_read;
	if ( avail <= 0 )
		return 0;
	
	int count = (avail - 1) / (period * 8) + 1;
	if ( !(regs [0] & loop_flag) && count > length_counter )
		count = length_counter;
	
	if ( last_read ) {
		*last_read = first_read + (count - 1) * (period * 8) + 1;
		assert( *last_read <= time );
		assert( count == count_reads( *last_read, NULL ) );
		assert( count - 1 == count_reads( *last_read - 1, NULL ) );
	}
	
	return count;
}

static const short dmc_period_table [2] [16] = {
	0x1ac, 0x17c, 0x154, 0x140, 0x11e, 0x0fe, 0x0e2, 0x0d6, // NTSC
	0x0be, 0x0a0, 0x08e, 0x080, 0x06a, 0x054, 0x048, 0x036,
	
	0x18e, 0x161, 0x13c, 0x129, 0x10a, 0x0ec, 0x0d2, 0x0c7, // PAL (totally untested)
	0x0b1, 0x095, 0x084, 0x077, 0x062, 0x04e, 0x043, 0x032  // to do: verify PAL periods
};

inline void Nes_Dmc::reload_sample()
{
	address = 0x4000 + regs [2] * 0x40;
	length_counter = regs [3] * 0x10 + 1;
}

static const unsigned char dac_table [128] = {
	 0,  0,  1,  2,  2,  3,  3,  4,  5,  5,  6,  7,  7,  8,  8,  9,
	10, 10, 11, 11, 12, 13, 13, 14, 14, 15, 15, 16, 17, 17, 18, 18,
	19, 19, 20, 20, 21, 21, 22, 22, 23, 23, 24, 24, 25, 25, 26, 26,
	27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 32, 32, 32, 33, 33, 34,
	34, 35, 35, 35, 36, 36, 37, 37, 38, 38, 38, 39, 39, 40, 40, 40,
	41, 41, 42, 42, 42, 43, 43, 44, 44, 44, 45, 45, 45, 46, 46, 47,
	47, 47, 48, 48, 48, 49, 49, 49, 50, 50, 50, 51, 51, 51, 52, 52,
	52, 53, 53, 53, 54, 54, 54, 55, 55, 55, 56, 56, 56, 57, 57, 57
};

void Nes_Dmc::write_register( int addr, int data )
{
	if ( addr == 0 ) {
		period = dmc_period_table [pal_mode] [data & 15];
		irq_enabled = (data & 0xc0) == 0x80; // enabled only if loop disabled
		irq_flag &= irq_enabled;
		recalc_irq();
	}
	else if ( addr == 1 )
	{
		if ( !nonlinear )
		{
			// adjust last_amp so that "pop" amplitude will be properly non-linear
			// with respect to change in dac
			int old_amp = dac_table [dac];
			dac = data & 0x7F;
			int diff = dac_table [dac] - old_amp;
			last_amp = dac - diff;
		}
		
		dac = data & 0x7F;
	}
}

void Nes_Dmc::start()
{
	reload_sample();
	fill_buffer();
	recalc_irq();
}

void Nes_Dmc::fill_buffer()
{
	if ( buf_empty && length_counter )
	{
		require( rom_reader ); // rom_reader must be set
		buf = rom_reader( rom_reader_data, 0x8000u + address );
		address = (address + 1) & 0x7FFF;
		buf_empty = false;
		if ( --length_counter == 0 )
		{
			if ( regs [0] & loop_flag ) {
				reload_sample();
			}
			else {
				apu->osc_enables &= ~0x10;
				irq_flag = irq_enabled;
				next_irq = Nes_Apu::no_irq;
				apu->irq_changed();
			}
		}
	}
}

void Nes_Dmc::run( cpu_time_t time, cpu_time_t end_time )
{
	if ( !output )
		return;
	
	int delta = update_amp( dac );
	if ( delta )
		synth.offset( time, delta, output );
	
	time += delay;
	if ( time < end_time )
	{
		int bits_remain = this->bits_remain;
		if ( silence && buf_empty )
		{
			int count = (end_time - time + period - 1) / period;
			bits_remain = (bits_remain - 1 + 8 - (count % 8)) % 8 + 1;
			time += count * period;
		}
		else
		{
			Blip_Buffer* const output = this->output;
			const int period = this->period;
			int bits = this->bits;
			int dac = this->dac;
			
			do
			{
				if ( !silence )
				{
					const int step = (bits & 1) * 4 - 2;
					bits >>= 1;
					if ( unsigned (dac + step) <= 0x7F ) {
						dac += step;
						synth.offset_inline( time, step, output );
					}
				}
				
				time += period;
				
				if ( --bits_remain == 0 )
				{
					bits_remain = 8;
					if ( buf_empty ) {
						silence = true;
					}
					else {
						silence = false;
						bits = buf;
						buf_empty = true;
						fill_buffer();
					}
				}
			}
			while ( time < end_time );
			
			this->dac = dac;
			this->last_amp = dac;
			this->bits = bits;
		}
		this->bits_remain = bits_remain;
	}
	delay = time - end_time;
}

// Nes_Noise

static const short noise_period_table [16] = {
	0x004, 0x008, 0x010, 0x020, 0x040, 0x060, 0x080, 0x0A0,
	0x0CA, 0x0FE, 0x17C, 0x1FC, 0x2FA, 0x3F8, 0x7F2, 0xFE4
};

void Nes_Noise::run( cpu_time_t time, cpu_time_t end_time )
{
	if ( !output )
		return;
	
	const int volume = this->volume();
	int amp = (noise & 1) ? volume : 0;
	int delta = update_amp( amp );
	if ( delta )
		synth.offset( time, delta, output );
	
	time += delay;
	if ( time < end_time )
	{
		const int mode_flag = 0x80;
		
		int period = noise_period_table [regs [2] & 15];
		if ( !volume )
		{
			// round to next multiple of period
			time += (end_time - time + period - 1) / period * period;
			
			// approximate noise cycling while muted, by shuffling up noise register
			// to do: precise muted noise cycling?
			if ( !(regs [2] & mode_flag) ) {
				int feedback = (noise << 13) ^ (noise << 14);
				noise = (feedback & 0x4000) | (noise >> 1);
			}
		}
		else
		{
			Blip_Buffer* const output = this->output;
			
			// using resampled time avoids conversion in synth.offset()
			Blip_Buffer::resampled_time_t rperiod = output->resampled_duration( period );
			Blip_Buffer::resampled_time_t rtime = output->resampled_time( time );
			
			int noise = this->noise;
			int delta = amp * 2 - volume;
			const int tap = (regs [2] & mode_flag ? 8 : 13);
			
			do {
				int feedback = (noise << tap) ^ (noise << 14);
				time += period;
				
				if ( (noise + 1) & 2 ) {
					// bits 0 and 1 of noise differ
					delta = -delta;
					synth.offset_resampled( rtime, delta, output );
				}
				
				rtime += rperiod;
				noise = (feedback & 0x4000) | (noise >> 1);
			}
			while ( time < end_time );
			
			last_amp = (delta + volume) >> 1;
			this->noise = noise;
		}
	}
	
	delay = time - end_time;
}

Blip_Buffer::Blip_Buffer()
{
	samples_per_sec = 44100;
	buffer_ = NULL;
	
	// try to cause assertion failure if buffer is used before these are set
	clocks_per_sec = 0;
	factor_ = ~0ul;
	offset_ = 0;
	buffer_size_ = 0;
	length_ = 0;
	
	bass_freq_ = 16;
}

void Blip_Buffer::clear( bool entire_buffer )
{
	long count = (entire_buffer ? buffer_size_ : samples_avail());
	offset_ = 0;
	reader_accum = 0;
	memset( buffer_, sample_offset & 0xFF, (count + widest_impulse_) * sizeof (buf_t_) );
}

blargg_err_t Blip_Buffer::sample_rate( long new_rate, int msec )
{
	unsigned new_size = (UINT_MAX >> BLIP_BUFFER_ACCURACY) + 1 - widest_impulse_ - 64;
	if ( msec != blip_default_length )
	{
		size_t s = (new_rate * (msec + 1) + 999) / 1000;
		if ( s < new_size )
			new_size = s;
		else
			require( false ); // requested buffer length exceeds limit
	}
	
	if ( buffer_size_ != new_size )
	{
		delete [] buffer_;
		buffer_ = NULL; // allow for exception in allocation below
		buffer_size_ = 0;
		offset_ = 0;
		
		buffer_ = BLARGG_NEW buf_t_ [new_size + widest_impulse_];
		BLARGG_CHECK_ALLOC( buffer_ );
	}
	
	buffer_size_ = new_size;
	length_ = new_size * 1000 / new_rate - 1;
	if ( msec )
		assert( length_ == msec ); // ensure length is same as that passed in
	
	samples_per_sec = new_rate;
	if ( clocks_per_sec )
		clock_rate( clocks_per_sec ); // recalculate factor
	
	bass_freq( bass_freq_ ); // recalculate shift
	
	clear();
	
	return blargg_success;
}

void Blip_Buffer::clock_rate( long cps )
{
	clocks_per_sec = cps;
	factor_ = (unsigned long) floor( (double) samples_per_sec / cps *
			(1L << BLIP_BUFFER_ACCURACY) + 0.5 );
	require( factor_ > 0 ); // clock_rate/sample_rate ratio is too large
}

Blip_Buffer::~Blip_Buffer()
{
	delete [] buffer_;
}

void Blip_Buffer::bass_freq( int freq )
{
	bass_freq_ = freq;
	if ( freq == 0 ) {
		bass_shift = 31; // 32 or greater invokes undefined behavior elsewhere
		return;
	}
	bass_shift = 1 + (int) floor( 1.442695041 * log( 0.124 * samples_per_sec / freq ) );
	if ( bass_shift < 0 )
		bass_shift = 0;
	if ( bass_shift > 24 )
		bass_shift = 24;
}

long Blip_Buffer::count_samples( blip_time_t t ) const {
	return (resampled_time( t ) >> BLIP_BUFFER_ACCURACY) - (offset_ >> BLIP_BUFFER_ACCURACY);
}

void Blip_Impulse_::init( blip_pair_t_* imps, int w, int r, int fb )
{
	fine_bits = fb;
	width = w;
	impulses = (imp_t*) imps;
	generate = true;
	volume_unit_ = -1.0;
	res = r;
	buf = NULL;
	
	impulse = &impulses [width * res * 2 * (fine_bits ? 2 : 1)];
	offset = 0;
}

const int impulse_bits = 15;
const long impulse_amp = 1L << impulse_bits;
const long impulse_offset = impulse_amp / 2;

void Blip_Impulse_::scale_impulse( int unit, imp_t* imp_in ) const
{
	long offset = ((long) unit << impulse_bits) - impulse_offset * unit +
			(1 << (impulse_bits - 1));
	imp_t* imp = imp_in;
	imp_t* fimp = impulse;
	for ( int n = res / 2 + 1; n--; )
	{
		int error = unit;
		for ( int nn = width; nn--; )
		{
			long a = ((long) *fimp++ * unit + offset) >> impulse_bits;
			error -= a - unit;
			*imp++ = (imp_t) a;
		}
		
		// add error to middle
		imp [-width / 2 - 1] += (imp_t) error;
	}
	
	if ( res > 2 ) {
		// second half is mirror-image
		const imp_t* rev = imp - width - 1;
		for ( int nn = (res / 2 - 1) * width - 1; nn--; )
			*imp++ = *--rev;
		*imp++ = (imp_t) unit;
	}
	
	// copy to odd offset
	*imp++ = (imp_t) unit;
	memcpy( imp, imp_in, (res * width - 1) * sizeof *imp );
}

const int max_res = 1 << blip_res_bits_;

void Blip_Impulse_::fine_volume_unit()
{
	// to do: find way of merging in-place without temporary buffer
	
	imp_t temp [max_res * 2 * Blip_Buffer::widest_impulse_];
	scale_impulse( (offset & 0xffff) << fine_bits, temp );
	imp_t* imp2 = impulses + res * 2 * width;
	scale_impulse( offset & 0xffff, imp2 );
	
	// merge impulses
	imp_t* imp = impulses;
	imp_t* src2 = temp;
	for ( int n = res / 2 * 2 * width; n--; ) {
		*imp++ = *imp2++;
		*imp++ = *imp2++;
		*imp++ = *src2++;
		*imp++ = *src2++;
	}
}

void Blip_Impulse_::volume_unit( double new_unit )
{
	if ( new_unit == volume_unit_ )
		return;
	
	if ( generate )
		treble_eq( blip_eq_t( -8.87, 8800, 44100 ) );
	
	volume_unit_ = new_unit;
	
	offset = 0x10001 * (unsigned long) floor( volume_unit_ * 0x10000 + 0.5 );
	
	if ( fine_bits )
		fine_volume_unit();
	else
		scale_impulse( offset & 0xffff, impulses );
}

static const double pi = 3.1415926535897932384626433832795029L;

void Blip_Impulse_::treble_eq( const blip_eq_t& new_eq )
{
	if ( !generate && new_eq.treble == eq.treble && new_eq.cutoff == eq.cutoff &&
			new_eq.sample_rate == eq.sample_rate )
		return; // already calculated with same parameters
	
	generate = false;
	eq = new_eq;
	
	double treble = pow( 10.0, 1.0 / 20 * eq.treble ); // dB (-6dB = 0.50)
	if ( treble < 0.000005 )
		treble = 0.000005;
	
	const double treble_freq = 22050.0; // treble level at 22 kHz harmonic
	const double sample_rate = eq.sample_rate;
	const double pt = treble_freq * 2 / sample_rate;
	double cutoff = eq.cutoff * 2 / sample_rate;
	if ( cutoff >= pt * 0.95 || cutoff >= 0.95 ) {
		cutoff = 0.5;
		treble = 1.0;
	}
	
	// DSF Synthesis (See T. Stilson & J. Smith (1996),
	// Alias-free digital synthesis of classic analog waveforms)
	
	// reduce adjacent impulse interference by using small part of wide impulse
	const double n_harm = 4096;
	const double rolloff = pow( treble, 1.0 / (n_harm * pt - n_harm * cutoff) );
	const double rescale = 1.0 / pow( rolloff, n_harm * cutoff );
	
	const double pow_a_n = rescale * pow( rolloff, n_harm );
	const double pow_a_nc = rescale * pow( rolloff, n_harm * cutoff );
	
	double total = 0.0;
	const double to_angle = pi / 2 / n_harm / max_res;
	
	float buf [max_res * (Blip_Buffer::widest_impulse_ - 2) / 2];
	const int size = max_res * (width - 2) / 2;
	for ( int i = size; i--; )
	{
		double angle = (i * 2 + 1) * to_angle;
		
		// equivalent
		//double y =     dsf( angle, n_harm * cutoff, 1.0 );
		//y -= rescale * dsf( angle, n_harm * cutoff, rolloff );
		//y += rescale * dsf( angle, n_harm,          rolloff );
		
		const double cos_angle = cos( angle );
		const double cos_nc_angle = cos( n_harm * cutoff * angle );
		const double cos_nc1_angle = cos( (n_harm * cutoff - 1.0) * angle );
		
		double b = 2.0 - 2.0 * cos_angle;
		double a = 1.0 - cos_angle - cos_nc_angle + cos_nc1_angle;
		
		double d = 1.0 + rolloff * (rolloff - 2.0 * cos_angle);
		double c = pow_a_n * rolloff * cos( (n_harm - 1.0) * angle ) -
				pow_a_n * cos( n_harm * angle ) -
				pow_a_nc * rolloff * cos_nc1_angle +
				pow_a_nc * cos_nc_angle;
		
		// optimization of a / b + c / d
		double y = (a * d + c * b) / (b * d);
		
		// fixed window which affects wider impulses more
		if ( width > 12 ) {
			double window = cos( n_harm / 1.25 / Blip_Buffer::widest_impulse_ * angle );
			y *= window * window;
		}
		
		total += (float) y;
		buf [i] = (float) y;
	}
	
	// integrate runs of length 'max_res'
	double factor = impulse_amp * 0.5 / total; // 0.5 accounts for other mirrored half
	imp_t* imp = impulse;
	const int step = max_res / res;
	int offset = res > 1 ? max_res : max_res / 2;
	for ( int n = res / 2 + 1; n--; offset -= step )
	{
		for ( int w = -width / 2; w < width / 2; w++ )
		{
			double sum = 0;
			for ( int i = max_res; i--; )
			{
				int index = w * max_res + offset + i;
				if ( index < 0 )
					index = -index - 1;
				if ( index < size )
					sum += buf [index];
			}
			*imp++ = (imp_t) floor( sum * factor + (impulse_offset + 0.5) );
		}
	}
	
	// rescale
	double unit = volume_unit_;
	if ( unit >= 0 ) {
		volume_unit_ = -1;
		volume_unit( unit );
	}
}

void Blip_Buffer::remove_samples( long count )
{
	require( buffer_ ); // sample rate must have been set
	
	if ( !count ) // optimization
		return;
	
	remove_silence( count );
	
	// Allows synthesis slightly past time passed to end_frame(), as long as it's
	// not more than an output sample.
	// to do: kind of hacky, could add run_until() which keeps track of extra synthesis
	int const copy_extra = 1;
	
	// copy remaining samples to beginning and clear old samples
	long remain = samples_avail() + widest_impulse_ + copy_extra;
	if ( count >= remain )
		memmove( buffer_, buffer_ + count, remain * sizeof (buf_t_) );
	else
		memcpy(  buffer_, buffer_ + count, remain * sizeof (buf_t_) );
	memset( buffer_ + remain, sample_offset & 0xFF, count * sizeof (buf_t_) );
}

long Blip_Buffer::read_samples( blip_sample_t* out, long max_samples, bool stereo )
{
	require( buffer_ ); // sample rate must have been set
	
	long count = samples_avail();
	if ( count > max_samples )
		count = max_samples;
	
	if ( !count )
		return 0; // optimization
	
	int sample_offset = this->sample_offset;
	int bass_shift = this->bass_shift;
	buf_t_* buf = buffer_;
	long accum = reader_accum;
	
	if ( !stereo ) {
		for ( long n = count; n--; ) {
			long s = accum >> accum_fract;
			accum -= accum >> bass_shift;
			accum += (long (*buf++) - sample_offset) << accum_fract;
			*out++ = (blip_sample_t) s;
			
			// clamp sample
			if ( (BOOST::int16_t) s != s )
				out [-1] = blip_sample_t (0x7FFF - (s >> 24));
		}
	}
	else {
		for ( long n = count; n--; ) {
			long s = accum >> accum_fract;
			accum -= accum >> bass_shift;
			accum += (long (*buf++) - sample_offset) << accum_fract;
			*out = (blip_sample_t) s;
			out += 2;
			
			// clamp sample
			if ( (BOOST::int16_t) s != s )
				out [-2] = blip_sample_t (0x7FFF - (s >> 24));
		}
	}
	
	reader_accum = accum;
	
	remove_samples( count );
	
	return count;
}

void Blip_Buffer::mix_samples( const blip_sample_t* in, long count )
{
	buf_t_* buf = &buffer_ [(offset_ >> BLIP_BUFFER_ACCURACY) + (widest_impulse_ / 2 - 1)];
	
	int prev = 0;
	while ( count-- ) {
		int s = *in++;
		*buf += s - prev;
		prev = s;
		++buf;
	}
	*buf -= *--in;
}

// Initial declarations
namespace APU {
	template <bool write> u8 access(int elapsed, u16 addr, u8 v = 0);
	void run_frame(int elapsed), reset(), init();
}
namespace CPU {
	enum IntType { NMI, RESET, IRQ, BRK }; // Interrupt type
	typedef u16 (*Mode)(void); // Addressing mode
	enum Flag { C, Z, I, D, V, N }; // Processor flags
	class Flags { bool f[6]; public:
		bool& operator[] (const int i) { return f[i]; }
		u8 get() { return f[C] | f[Z] << 1 | f[I] << 2 | f[D] << 3 | 1 << 5 | f[V] << 6 | f[N] << 7; }
		void set(u8 p) { f[C] = NTH_BIT(p, 0), f[Z] = NTH_BIT(p, 1), f[I] = NTH_BIT(p, 2), f[D] = NTH_BIT(p, 3), f[V] = NTH_BIT(p, 6), f[N] = NTH_BIT(p, 7); }
	};
	int dmc_read(void*, cpu_addr_t addr);
	void set_nmi(bool v = true), set_irq(bool v = true), power(), run_frame();
	// Save states
	struct save_state {	u8 A, X, Y, S, ram[0x800]; u16 PC; CPU::Flags P; bool nmi, irq; };
	void save(), load();
}
namespace PPU {
	enum Scanline  { VISIBLE, POST, NMI, PRE };
	enum Mirroring { VERTICAL, HORIZONTAL, ONE_SCREEN_HI, ONE_SCREEN_LO, FOUR_SCREEN };
	struct Sprite { // Sprite buffer
		// Index in OAM, X position, Y position, Tile index, Attributes, Tile data (low), Tile data (high)
		u8 id, x, y, tile, attr, dataL, dataH;
	};
	union Ctrl { // PPUCTRL ($2000) register
		struct {
			// Nametable, Address increment, Sprite pattern table, BG pattern table, Sprite size, PPU master/slave, Enable NMI 
			unsigned nt : 2, incr : 1, sprTbl : 1, bgTbl : 1, sprSz : 1, slave : 1, nmi : 1;
		};
		u8 r;
	};
	union Mask { // PPUMASK ($2001) register
		struct {
			// Grayscale; Show background in leftmost 8 pixels, sprite in leftmost 8 pixels, background, sprites; Intensify reds, greens, blues
			unsigned gray : 1, bgLeft : 1, sprLeft : 1, bg : 1, spr : 1, red : 1, green : 1, blue : 1;
		};
		u8 r;
	};
	union Status { // PPUSTATUS ($2002) register
		struct {
			// Not significant, Sprite overflow, Sprite 0 Hit, In VBlank
			unsigned bus : 5, sprOvf : 1, sprHit : 1, vBlank : 1;
		};
		u8 r;
	};
	union Addr { // Loopy's VRAM address
		struct {
			// Coarse X, Coarse Y, Nametable, Fine Y
			unsigned cX : 5, cY : 5, nt : 2, fY : 3;
		};
		struct { unsigned l : 8, h : 7; };
		unsigned addr : 14, r : 15;
	};
	template <bool write> u8 access(u16 index, u8 v = 0);
	void set_mirroring(Mirroring mode), step(), reset();
	// Save states
	struct save_state { u8 ciRam[0x800], cgRam[0x20], oamMem[0x100]; Sprite oam[8], secOam[8]; u32 pixels[256 * 240]; Ctrl ctrl; Mask mask; Status status; };
	void save(), load();
}
namespace Cartridge {
	template <bool wr> u8 access(u16 addr, u8 v = 0); template <bool wr> u8 chr_access(u16 addr, u8 v = 0);
	void signal_scanline(), load(const char *fileName);
}
namespace Joypad {
	u8 read_state(int n); void write_strobe(bool v);
}
namespace GUI {
	const unsigned WIDTH  = 256, HEIGHT = 240; // Screen size
	int query_button();
	void init(), run(), save(), load(), new_frame(u32* pixels), new_samples(const blip_sample_t* samples, size_t count);;
	u8 get_joypad_state(int n);
	SDL_Scancode query_key();
}

// Mappers
class Mapper {
	u8* rom; bool chrRam = false;
protected:
	u8 *prg, *chr, *prgRam;
	u32 prgSize, chrSize, prgRamSize, prgMap[4], chrMap[8];
	template <int pageKBs> void map_prg(int slot, int bank);
	template <int pageKBs> void map_chr(int slot, int bank);
public:
	Mapper(u8* rom);
	~Mapper();
	u8 read(u16 addr), chr_read(u16 addr);
	virtual u8 write(u16 addr, u8 v) { return v; }
	virtual u8 chr_write(u16 addr, u8 v) { return v; }
	virtual void signal_scanline() {}
};
Mapper::Mapper(u8* rom) : rom(rom) {
	// Read infos from header
	prgSize = rom[4] * 0x4000, chrSize = rom[5] * 0x2000;
	prgRamSize = rom[8] ? rom[8] * 0x2000 : 0x2000;
	set_mirroring((rom[6] & 1) ? PPU::VERTICAL : PPU::HORIZONTAL);
	this->prg = rom + 16, this->prgRam = new u8[prgRamSize];
	if (chrSize) this->chr = rom + 16 + prgSize; // CHR ROM
	else { chrRam = true, chrSize = 0x2000, this->chr = new u8[chrSize]; } // CHR RAM
}
Mapper::~Mapper() { delete rom; delete prgRam; if (chrRam) delete chr; }
// Access to memory
u8 Mapper::read(u16 addr) { return addr >= 0x8000 ? prg[prgMap[(addr - 0x8000) / 0x2000] + ((addr - 0x8000) % 0x2000)] : prgRam[addr - 0x6000]; }
u8 Mapper::chr_read(u16 addr) { return chr[chrMap[addr / 0x400] + (addr % 0x400)]; }
template <int pageKBs> void Mapper::map_prg(int slot, int bank) { // PRG mapping functions
	if (bank < 0) bank = (prgSize / (0x400*pageKBs)) + bank;
	for (int i = 0; i < (pageKBs/8); i++) prgMap[(pageKBs/8) * slot + i] = (pageKBs*0x400*bank + 0x2000*i) % prgSize;
}
template void Mapper::map_prg<32>(int, int); template void Mapper::map_prg<16>(int, int); template void Mapper::map_prg<8> (int, int);
template <int pageKBs> void Mapper::map_chr(int slot, int bank) { // CHR mapping functions
	for (int i = 0; i < pageKBs; i++) chrMap[pageKBs*slot + i] = (pageKBs*0x400*bank + 0x400*i) % chrSize;
}
template void Mapper::map_chr<8>(int, int); template void Mapper::map_chr<4>(int, int); template void Mapper::map_chr<2>(int, int); template void Mapper::map_chr<1>(int, int);
class Mapper0 : public Mapper {
	public: Mapper0(u8* rom) : Mapper(rom) { map_prg<32>(0, 0); map_chr<8> (0, 0); }
};
class Mapper1 : public Mapper {
	int writeN; u8 tmpReg, regs[4];
	void apply() { // Apply the registers state 
		if (regs[0] & 0b1000) { // 16KB PRG
			// 0x8000 swappable, 0xC000 fixed to bank 0x0F
			if (regs[0] & 0b100) { map_prg<16>(0, regs[3] & 0xF); map_prg<16>(1, 0xF); }
			// 0x8000 fixed to bank 0x00, 0xC000 swappable
			else { map_prg<16>(0, 0); map_prg<16>(1, regs[3] & 0xF); }
		}
		else map_prg<32>(0, (regs[3] & 0xF) >> 1); // 32KB PRG
		if (regs[0] & 0b10000) { map_chr<4>(0, regs[1]); map_chr<4>(1, regs[2]); } // 4KB CHR
		else map_chr<8>(0, regs[1] >> 1); // 8KB CHR
		switch (regs[0] & 0b11) { // Set mirroring
			case 0: set_mirroring(PPU::ONE_SCREEN_LO); break;
			case 1: set_mirroring(PPU::ONE_SCREEN_HI); break;
        	case 2: set_mirroring(PPU::VERTICAL);      break;
        	case 3: set_mirroring(PPU::HORIZONTAL);    break;
		}
	}
public:
 	Mapper1(u8* rom) : Mapper(rom) {
		regs[0] = 0x0C;
		writeN = tmpReg = regs[1] = regs[2] = regs[3] = 0;
		apply();
	}
	u8 write(u16 addr, u8 v) {
		if (addr < 0x8000) prgRam[addr - 0x6000] = v; // PRG RAM write
		else if (addr & 0x8000) { // Mapper register write
			if (v & 0x80) { writeN = tmpReg = 0; regs[0] |= 0x0C; apply(); } // Reset
			else {
				tmpReg = ((v & 1) << 4) | (tmpReg >> 1); // Write a bit into the temporary register
				if (++writeN == 5) { regs[(addr >> 13) & 0b11] = tmpReg; writeN = tmpReg = 0; apply(); } // Finished writing all the bits
			}
		}
		return v;
	}
	u8 chr_write(u16 addr, u8 v) { return chr[addr] = v; }
};
class Mapper2 : public Mapper {
	u8 regs[1]; bool vertical_mirroring;
	void apply() { // Apply the registers state
		// 16 kb PRG ROM Banks: 0x8000 - 0xBFFF swappable, 0xC000 - 0xFFFF fixed
		map_prg<16>(0, regs[0] & 0xF); map_prg<16>(1, 0xF);
		map_chr<8>(0, 0); // 8k of CHR
		set_mirroring(vertical_mirroring?PPU::VERTICAL:PPU::HORIZONTAL); // mirroring is based on the header (soldered)
	}
public:
	Mapper2(u8* rom) : Mapper(rom) {
		regs[0] = 0;
		vertical_mirroring = rom[6] & 0x01;
		apply();
	}
	u8 write(u16 addr, u8 v) {
		if (addr & 0x8000) { regs[0] = v; apply(); } // bank switching
		return v;
	}
	u8 chr_write(u16 addr, u8 v) { return chr[addr] = v; }
};
class Mapper3 : public Mapper {
	u8 regs[1]; bool vertical_mirroring, PRG_size_16k;
	void apply() { // Apply the registers state
		if (PRG_size_16k) {	map_prg<16>(0,0); map_prg<16>(1,0); } // mirror the bottom on the top: 0x8000 - 0xBFFF == 0xC000 - 0xFFFF
		else { map_prg<16>(0,0); map_prg<16>(1,1); } // no mirroring
		map_chr<8>(0, regs[0] & 0b11); // 8k bankswitched CHR		
		set_mirroring(vertical_mirroring?PPU::VERTICAL:PPU::HORIZONTAL); // mirroring is based on the header (soldered)
	}
public:
	Mapper3(u8* rom) : Mapper(rom) {
		PRG_size_16k = rom[4] == 1;
		vertical_mirroring = rom[6] & 0x01;
		regs[0] = 0;
		apply();
	}
	u8 write(u16 addr, u8 v) {
		if (addr & 0x8000) { regs[0] = v; apply(); } // chr bank switching
		return v;
	}
	u8 chr_write(u16 addr, u8 v) { return chr[addr] = v; }
};
class Mapper4 : public Mapper {
	u8 reg8000, regs[8], irqPeriod, irqCounter; bool horizMirroring, irqEnabled;
	void apply() { // Apply the registers state
		map_prg<8>(1, regs[7]);
		if (!(reg8000 & (1 << 6))) { map_prg<8>(0, regs[6]); map_prg<8>(2, -2); } // PRG Mode 0
		else { map_prg<8>(0, -2); map_prg<8>(2, regs[6]); } // PRG Mode 1
		if (!(reg8000 & (1 << 7))) { // CHR Mode 0
			map_chr<2>(0, regs[0] >> 1); map_chr<2>(1, regs[1] >> 1);
			for (int i = 0; i < 4; i++) map_chr<1>(4 + i, regs[2 + i]);
		}
		else { // CHR Mode 1
			for (int i = 0; i < 4; i++) map_chr<1>(i, regs[2 + i]);
			map_chr<2>(2, regs[0] >> 1); map_chr<2>(3, regs[1] >> 1);
		}
		set_mirroring(horizMirroring ? PPU::HORIZONTAL : PPU::VERTICAL);
	}
public:
	Mapper4(u8* rom) : Mapper(rom) {
		for (int i = 0; i < 8; i++) regs[i] = 0;
		horizMirroring = true; irqEnabled = false; irqPeriod = irqCounter = 0;
		map_prg<8>(3, -1);
		apply();
	}
	u8 write(u16 addr, u8 v) {
		if (addr < 0x8000) prgRam[addr - 0x6000] = v;
		else if (addr & 0x8000) {
			switch (addr & 0xE001) {
				case 0x8000:  reg8000 = v;                      break;
				case 0x8001:  regs[reg8000 & 0b111] = v;        break;
				case 0xA000:  horizMirroring = v & 1;           break;
				case 0xC000:  irqPeriod = v;                    break;
				case 0xC001:  irqCounter = 0;                   break;
				case 0xE000:  CPU::set_irq(irqEnabled = false); break;
				case 0xE001:  irqEnabled = true;                break;
			}
			apply();
		}
		return v;
	}
	u8 chr_write(u16 addr, u8 v) { return chr[addr] = v; }
	void signal_scanline() {
		irqCounter == 0 ? irqCounter = irqPeriod : --irqCounter;
		if (irqEnabled and irqCounter == 0)	CPU::set_irq();
	}
};
class Mapper7 : public Mapper {
    u8 regs[1];
	
    void apply() { // Apply the registers state
    	map_prg<32>(0, regs[0] & 0b00001111); // 32 kb PRG ROM Banks, 0x8000 - 0xFFFF swappable
    	map_chr<8>(0, 0); // 8k of CHR (ram)
    	set_mirroring((regs[0] & 0b00010000) ? PPU::ONE_SCREEN_HI : PPU::ONE_SCREEN_LO); // Mirroring based on bit 5
	}
public:
    Mapper7(u8* rom) : Mapper(rom) { regs[0] = 0; apply(); }
	u8 write(u16 addr, u8 v) {
    	if (addr & 0x8000) { regs[0] = v; apply(); } // Bank switching
    	return v;
	}
	u8 chr_write(u16 addr, u8 v) { return chr[addr] = v; }
};

// Actual code
namespace APU {
	Nes_Apu apu; Blip_Buffer buf;
	const int OUT_SIZE = 4096;
	blip_sample_t outBuf[OUT_SIZE];
	void init() {
		buf.sample_rate(96000); buf.clock_rate(1789773);
		apu.output(&buf); apu.dmc_reader(CPU::dmc_read);
	}
	void reset() { apu.reset(); buf.clear(); }
	template <bool write> u8 access(int elapsed, u16 addr, u8 v) {
		if (write) apu.write_register(elapsed, addr, v);
		else if (addr == apu.status_addr) v = apu.read_status(elapsed);
		return v;
	}
	template u8 access<0>(int, u16, u8); template u8 access<1>(int, u16, u8);
	void run_frame(int elapsed) {
		apu.end_frame(elapsed); buf.end_frame(elapsed);
		if (buf.samples_avail() >= OUT_SIZE) GUI::new_samples(outBuf, buf.read_samples(outBuf, OUT_SIZE));
	}
}
namespace CPU {
	u8 A, X, Y, S, ram[0x800]; u16 PC; Flags P; bool nmi, irq; // CPU state
	save_state ss; // Save state
	// Remaining clocks to end frame
	const int TOTAL_CYCLES = 29781; int remainingCycles;
	inline int elapsed() { return TOTAL_CYCLES - remainingCycles; }
	// Cycle emulation
	#define T tick()
	inline void tick() { PPU::step(); PPU::step(); PPU::step(); remainingCycles--; }
	// Flags updating
	inline void upd_cv(u8 x, u8 y, s16 r) { P[C] = (r>0xFF); P[V] = ~(x^y) & (x^r) & 0x80; }
	inline void upd_nz(u8 x) { P[N] = x & 0x80; P[Z] = (x == 0); }
	inline bool cross(u16 a, u8 i) { return ((a+i) & 0xFF00) != ((a & 0xFF00)); } // Does adding I to A cross a page?
	// Memory access
	void dma_oam(u8 bank);
	template<bool wr> inline u8 access(u16 addr, u8 v = 0) {
		u8* r;
		switch (addr) {
			case 0x0000 ... 0x1FFF:  r = &ram[addr % 0x800]; if (wr) *r = v; return *r; // RAM
			case 0x2000 ... 0x3FFF:  return PPU::access<wr>(addr % 8, v);               // PPU
			case 0x4000 ... 0x4013:                                                     // APU 
			case            0x4015:          return APU::access<wr>(elapsed(), addr, v);
			case            0x4017:  if (wr) return APU::access<wr>(elapsed(), addr, v);
									 else return Joypad::read_state(1);              // Joypad 1
			case            0x4014:  if (wr) dma_oam(v); break;                      // OAM DMA
			case            0x4016:  if (wr) { Joypad::write_strobe(v & 1); break; } // Joypad strobe
									 else return Joypad::read_state(0);              // Joypad 0
			case 0x4018 ... 0xFFFF:  return Cartridge::access<wr>(addr, v);          // Cartridge
		}
		return 0;
	}
	inline u8  wr(u16 a, u8 v)      { T; return access<1>(a, v);   }
	inline u8  rd(u16 a)            { T; return access<0>(a);      }
	inline u16 rd16_d(u16 a, u16 b) { return rd(a) | (rd(b) << 8); } // Read from A and B and merge
	inline u16 rd16(u16 a)          { return rd16_d(a, a+1);       }
	inline u8  push(u8 v)           { return wr(0x100 + (S--), v); }
	inline u8  pop()                { return rd(0x100 + (++S));    }
	void dma_oam(u8 bank) { for (int i = 0; i < 256; i++)  wr(0x2014, rd(bank*0x100 + i)); }
	// Addressing modes
	inline u16 imm()   { return PC++;                                       }
	inline u16 imm16() { PC += 2; return PC - 2;                            }
	inline u16 abs()   { return rd16(imm16());                              }
	inline u16 _abx()  { T; return abs() + X;                               } // Exception
	inline u16 abx()   { u16 a = abs(); if (cross(a, X)) T; return a + X;   }
	inline u16 aby()   { u16 a = abs(); if (cross(a, Y)) T; return a + Y;   }
	inline u16 zp()    { return rd(imm());                                  }
	inline u16 zpx()   { T; return (zp() + X) % 0x100;                      }
	inline u16 zpy()   { T; return (zp() + Y) % 0x100;                      }
	inline u16 izx()   { u8 i = zpx(); return rd16_d(i, (i+1) % 0x100);     }
	inline u16 _izy()  { u8 i = zp();  return rd16_d(i, (i+1) % 0x100) + Y; } // Exception
	inline u16 izy()   { u16 a = _izy(); if (cross(a-Y, Y)) T; return a;    }
	// STx
	template<u8& r, Mode m> void st()        {    wr(   m()    , r); }
	template<>              void st<A,izy>() { T; wr(_izy()    , A); } // Exceptions
	template<>              void st<A,abx>() { T; wr( abs() + X, A); } // ...
	template<>              void st<A,aby>() { T; wr( abs() + Y, A); } // ...
	#define G  u16 a = m(); u8 p = rd(a) // Fetch parameter
	template<u8& r, Mode m> void ld()  { G; upd_nz(r = p);                  } // LDx
	template<u8& r, Mode m> void cmp() { G; upd_nz(r - p); P[C] = (r >= p); } // CMP, CPx
	// Arithmetic and bitwise
	template<Mode m> void ADC() { G       ; s16 r = A + p + P[C]; upd_cv(A, p, r); upd_nz(A = r); }
	template<Mode m> void SBC() { G ^ 0xFF; s16 r = A + p + P[C]; upd_cv(A, p, r); upd_nz(A = r); }
	template<Mode m> void BIT() { G; P[Z] = !(A & p); P[N] = p & 0x80; P[V] = p & 0x40; }
	template<Mode m> void AND() { G; upd_nz(A &= p); }
	template<Mode m> void EOR() { G; upd_nz(A ^= p); }
	template<Mode m> void ORA() { G; upd_nz(A |= p); }
	// Read-Modify-Write
	template<Mode m> void ASL() { G; P[C] = p & 0x80; T; upd_nz(wr(a, p << 1)); }
	template<Mode m> void LSR() { G; P[C] = p & 0x01; T; upd_nz(wr(a, p >> 1)); }
	template<Mode m> void ROL() { G; u8 c = P[C]     ; P[C] = p & 0x80; T; upd_nz(wr(a, (p << 1) | c) ); }
	template<Mode m> void ROR() { G; u8 c = P[C] << 7; P[C] = p & 0x01; T; upd_nz(wr(a, c | (p >> 1)) ); }
	template<Mode m> void DEC() { G; T; upd_nz(wr(a, --p)); }
	template<Mode m> void INC() { G; T; upd_nz(wr(a, ++p)); }
	#undef G
	// DEx, INx
	template<u8& r> void dec() { upd_nz(--r); T; }
	template<u8& r> void inc() { upd_nz(++r); T; }
	// Bit shifting on the accumulator
	void ASL_A() { P[C] = A & 0x80; upd_nz(A <<= 1); T; }
	void LSR_A() { P[C] = A & 0x01; upd_nz(A >>= 1); T; }
	void ROL_A() { u8 c = P[C]     ; P[C] = A & 0x80; upd_nz(A = ((A << 1) | c) ); T; }
	void ROR_A() { u8 c = P[C] << 7; P[C] = A & 0x01; upd_nz(A = (c | (A >> 1)) ); T; }
	// Txx (move values between registers)
	template<u8& s, u8& d> void tr()      { upd_nz(d = s); T; }
	template<>             void tr<X,S>() { S = X;         T; } // TSX, exception
	// Stack operations
	void PLP() { T; T; P.set(pop()); }
	void PHP() { T; push(P.get() | (1 << 4)); } // B flag set.
	void PLA() { T; T; A = pop(); upd_nz(A);  }
	void PHA() { T; push(A); }
	template<Flag f, bool v> void br() { // Flow control (branches, jumps)
		s8 j = rd(imm()); 
		if (P[f] == v) {
			if (cross(PC, j)) T;
			T; PC += j; 
		}
	}
	void JMP_IND() { u16 i = rd16(imm16()); PC = rd16_d(i, (i&0xFF00) | ((i+1) % 0x100)); }
	void JMP()     { PC = rd16(imm16()); }
	void JSR()     { u16 t = PC+1; T; push(t >> 8); push(t); PC = rd16(imm16()); }
	// Return instructions
	void RTS() { T; T;  PC = (pop() | (pop() << 8)) + 1; T; }
	void RTI() { PLP(); PC =  pop() | (pop() << 8);         }
	template<Flag f, bool v> void flag() { P[f] = v; T; } // Clear and set flags.
	template<IntType t> void INT() {
		T; if (t != BRK) T; // BRK already performed the fetch
		if (t != RESET) { // Writes on stack are inhibited on RESET
			push(PC >> 8); push(PC & 0xFF);
			push(P.get() | ((t == BRK) << 4));  // Set B if BRK
		}
		else { S -= 3; T; T; T; }
		P[I] = true; 			// NMI    Reset    IRQ     BRK  
		constexpr u16 vect[] = { 0xFFFA, 0xFFFC, 0xFFFE, 0xFFFE };
		PC = rd16(vect[t]);
		if (t == NMI) nmi = false;
	}
	void NOP() { T; }
	void exec() { // Execute a CPU instruction
		switch (rd(PC++)) { // Fetch the opcode and select the right function to emulate the instruction:
			case 0x00: return INT<BRK>()  ; case 0x01: return ORA<izx>()  ; case 0x05: return ORA<zp>()   ; case 0x06: return ASL<zp>()   ;
			case 0x08: return PHP()       ; case 0x09: return ORA<imm>()  ;	case 0x0A: return ASL_A()     ; case 0x0D: return ORA<abs>()  ;
			case 0x0E: return ASL<abs>()  ; case 0x10: return br<N,0>()   ;	case 0x11: return ORA<izy>()  ; case 0x15: return ORA<zpx>()  ;
			case 0x16: return ASL<zpx>()  ; case 0x18: return flag<C,0>() ;	case 0x19: return ORA<aby>()  ; case 0x1D: return ORA<abx>()  ;
			case 0x1E: return ASL<_abx>() ; case 0x20: return JSR()       ;	case 0x21: return AND<izx>()  ; case 0x24: return BIT<zp>()   ;
			case 0x25: return AND<zp>()   ; case 0x26: return ROL<zp>()   ;	case 0x28: return PLP()       ; case 0x29: return AND<imm>()  ;
			case 0x2A: return ROL_A()     ; case 0x2C: return BIT<abs>()  ;	case 0x2D: return AND<abs>()  ; case 0x2E: return ROL<abs>()  ;
			case 0x30: return br<N,1>()   ; case 0x31: return AND<izy>()  ;	case 0x35: return AND<zpx>()  ; case 0x36: return ROL<zpx>()  ;
			case 0x38: return flag<C,1>() ; case 0x39: return AND<aby>()  ;	case 0x3D: return AND<abx>()  ; case 0x3E: return ROL<_abx>() ;
			case 0x40: return RTI()       ; case 0x41: return EOR<izx>()  ;	case 0x45: return EOR<zp>()   ; case 0x46: return LSR<zp>()   ;
			case 0x48: return PHA()       ; case 0x49: return EOR<imm>()  ;	case 0x4A: return LSR_A()     ; case 0x4C: return JMP()       ;
			case 0x4D: return EOR<abs>()  ; case 0x4E: return LSR<abs>()  ; case 0x50: return br<V,0>()   ; case 0x51: return EOR<izy>()  ;
			case 0x55: return EOR<zpx>()  ; case 0x56: return LSR<zpx>()  ; case 0x58: return flag<I,0>() ; case 0x59: return EOR<aby>()  ;
			case 0x5D: return EOR<abx>()  ; case 0x5E: return LSR<_abx>() ;	case 0x60: return RTS()       ; case 0x61: return ADC<izx>()  ;
			case 0x65: return ADC<zp>()   ; case 0x66: return ROR<zp>()   ;	case 0x68: return PLA()       ; case 0x69: return ADC<imm>()  ;
			case 0x6A: return ROR_A()     ; case 0x6C: return JMP_IND()   ;	case 0x6D: return ADC<abs>()  ; case 0x6E: return ROR<abs>()  ;
			case 0x70: return br<V,1>()   ; case 0x71: return ADC<izy>()  ;	case 0x75: return ADC<zpx>()  ; case 0x76: return ROR<zpx>()  ;
			case 0x78: return flag<I,1>() ; case 0x79: return ADC<aby>()  ;	case 0x7D: return ADC<abx>()  ; case 0x7E: return ROR<_abx>() ;
			case 0x81: return st<A,izx>() ; case 0x84: return st<Y,zp>()  ;	case 0x85: return st<A,zp>()  ; case 0x86: return st<X,zp>()  ;
			case 0x88: return dec<Y>()    ; case 0x8A: return tr<X,A>()   ;	case 0x8C: return st<Y,abs>() ; case 0x8D: return st<A,abs>() ;
			case 0x8E: return st<X,abs>() ; case 0x90: return br<C,0>()   ;	case 0x91: return st<A,izy>() ; case 0x94: return st<Y,zpx>() ;
			case 0x95: return st<A,zpx>() ; case 0x96: return st<X,zpy>() ;	case 0x98: return tr<Y,A>()   ; case 0x99: return st<A,aby>() ;
			case 0x9A: return tr<X,S>()   ; case 0x9D: return st<A,abx>() ;	case 0xA0: return ld<Y,imm>() ; case 0xA1: return ld<A,izx>() ;
			case 0xA2: return ld<X,imm>() ; case 0xA4: return ld<Y,zp>()  ;	case 0xA5: return ld<A,zp>()  ; case 0xA6: return ld<X,zp>()  ;
			case 0xA8: return tr<A,Y>()   ; case 0xA9: return ld<A,imm>() ;	case 0xAA: return tr<A,X>()   ; case 0xAC: return ld<Y,abs>() ;
			case 0xAD: return ld<A,abs>() ; case 0xAE: return ld<X,abs>() ;	case 0xB0: return br<C,1>()   ; case 0xB1: return ld<A,izy>() ;
			case 0xB4: return ld<Y,zpx>() ; case 0xB5: return ld<A,zpx>() ;	case 0xB6: return ld<X,zpy>() ; case 0xB8: return flag<V,0>() ;
			case 0xB9: return ld<A,aby>() ; case 0xBA: return tr<S,X>()   ;	case 0xBC: return ld<Y,abx>() ; case 0xBD: return ld<A,abx>() ;
			case 0xBE: return ld<X,aby>() ; case 0xC0: return cmp<Y,imm>();	case 0xC1: return cmp<A,izx>(); case 0xC4: return cmp<Y,zp>() ;
			case 0xC5: return cmp<A,zp>() ; case 0xC6: return DEC<zp>()   ;	case 0xC8: return inc<Y>()    ; case 0xC9: return cmp<A,imm>();
			case 0xCA: return dec<X>()    ; case 0xCC: return cmp<Y,abs>();	case 0xCD: return cmp<A,abs>(); case 0xCE: return DEC<abs>()  ;
			case 0xD0: return br<Z,0>()   ; case 0xD1: return cmp<A,izy>();	case 0xD5: return cmp<A,zpx>(); case 0xD6: return DEC<zpx>()  ;
			case 0xD8: return flag<D,0>() ; case 0xD9: return cmp<A,aby>();	case 0xDD: return cmp<A,abx>(); case 0xDE: return DEC<_abx>() ;
			case 0xE0: return cmp<X,imm>(); case 0xE1: return SBC<izx>()  ;	case 0xE4: return cmp<X,zp>() ; case 0xE5: return SBC<zp>()   ;
			case 0xE6: return INC<zp>()   ; case 0xE8: return inc<X>()    ;	case 0xE9: return SBC<imm>()  ; case 0xEA: return NOP()       ;
			case 0xEC: return cmp<X,abs>(); case 0xED: return SBC<abs>()  ;	case 0xEE: return INC<abs>()  ; case 0xF0: return br<Z,1>()   ;
			case 0xF1: return SBC<izy>()  ; case 0xF5: return SBC<zpx>()  ;	case 0xF6: return INC<zpx>()  ; case 0xF8: return flag<D,1>() ;
			case 0xF9: return SBC<aby>()  ; case 0xFD: return SBC<abx>()  ;	case 0xFE: return INC<_abx>() ; default:
			std::cout << "Invalid Opcode! PC: " << PC << " Opcode: 0x" << std::hex << (int)(rd(PC - 1)) << "\n";
			return NOP();
		}
	}
	void set_nmi(bool v) { nmi = v; }
	void set_irq(bool v) { irq = v; }
	int dmc_read(void*, cpu_addr_t addr) { return access<0>(addr); }
	void power() { // Turn on the CPU
		remainingCycles = 0;
		P.set(0x04);
		A = X = Y = S = 0x00;
		memset(ram, 0xFF, sizeof(ram));
		nmi = irq = false;
		INT<RESET>();
	}
	void run_frame() { // Run the CPU for roughly a frame
		remainingCycles += TOTAL_CYCLES;
		while (remainingCycles > 0) {
			if (nmi) INT<NMI>();
			else if (irq and !P[I]) INT<IRQ>();
			exec();
		}
		APU::run_frame(elapsed());
	}
	void save() { // Save state
		ss.A = A, ss.X = X, ss.Y = Y, ss.S = S, ss.PC = PC, ss.P = P, ss.nmi = nmi, ss.irq = irq;
		memcpy(ss.ram, ram, sizeof ram);
	}
	void load() { // Load state
		A = ss.A, X = ss.X, Y = ss.Y, S = ss.S, PC = ss.PC, P = ss.P, nmi = ss.nmi, irq = ss.irq;
		memcpy(ram, ss.ram, sizeof ram);
	}	
}
namespace PPU {
	u32 nesRgb[] = { 0x7C7C7C, 0x0000FC, 0x0000BC, 0x4428BC, 0x940084, 0xA80020, 0xA81000, 0x881400,
					 0x503000, 0x007800, 0x006800, 0x005800, 0x004058, 0x000000, 0x000000, 0x000000,
					 0xBCBCBC, 0x0078F8, 0x0058F8, 0x6844FC, 0xD800CC, 0xE40058, 0xF83800, 0xE45C10,
					 0xAC7C00, 0x00B800, 0x00A800, 0x00A844, 0x008888, 0x000000, 0x000000, 0x000000,
					 0xF8F8F8, 0x3CBCFC, 0x6888FC, 0x9878F8, 0xF878F8, 0xF85898, 0xF87858, 0xFCA044,
					 0xF8B800, 0xB8F818, 0x58D854, 0x58F898, 0x00E8D8, 0x787878, 0x000000, 0x000000,
					 0xFCFCFC, 0xA4E4FC, 0xB8B8F8, 0xD8B8F8, 0xF8B8F8, 0xF8A4C0, 0xF0D0B0, 0xFCE0A8,
					 0xF8D878, 0xD8F878, 0xB8F8B8, 0xB8F8D8, 0x00FCFC, 0xF8D8F8, 0x000000, 0x000000 };
	Mirroring mirroring; // Mirroring mode
	u8 ciRam[0x800], cgRam[0x20], oamMem[0x100]; // VRAM for nametables, palettes, sprite properties
	Sprite oam[8], secOam[8]; // Sprite buffers
	u32 pixels[256 * 240]; // Video buffer
	Addr vAddr, tAddr; // Loopy V, T
	u8 fX, oamAddr; // Fine X, OAM address
	Ctrl ctrl;     // PPUCTRL   ($2000) register
	Mask mask;     // PPUMASK   ($2001) register
	Status status; // PPUSTATUS ($2002) register
	u8 nt, at, bgL, bgH, atShiftL, atShiftH; u16 bgShiftL, bgShiftH; // Background latches, shift registers
	bool atLatchL, atLatchH;
	int scanline, dot; bool frameOdd; // Rendering counters
	save_state ss;
	inline bool rendering() { return mask.bg || mask.spr; }
	inline int spr_height() { return ctrl.sprSz ? 16 : 8; }
	u16 nt_mirror(u16 addr) { // Get CIRAM address according to mirroring
		switch (mirroring) {
			case VERTICAL:   return addr % 0x800;
			case HORIZONTAL: return ((addr / 2) & 0x400) + (addr % 0x400);
			case ONE_SCREEN_LO:
        	case ONE_SCREEN_HI: return ((addr & 0x3ff) + ((mirroring == ONE_SCREEN_HI) ? 0x400 : 0x0)) - 0x2000;
			default:         return addr - 0x2000;
		}
	}
	void set_mirroring(Mirroring mode) { mirroring = mode; }
	u8 rd(u16 addr) { // Access PPU memory
		switch (addr) {
			case 0x0000 ... 0x1FFF: return Cartridge::chr_access<0>(addr);  // CHR-ROM/RAM
			case 0x2000 ... 0x3EFF: return ciRam[nt_mirror(addr)];          // Nametables
			case 0x3F00 ... 0x3FFF: // Palettes
				if ((addr & 0x13) == 0x10) addr &= ~0x10;
				return cgRam[addr & 0x1F] & (mask.gray ? 0x30 : 0xFF);
			default: return 0;
		}
	}
	void wr(u16 addr, u8 v) {
		switch (addr) {
			case 0x0000 ... 0x1FFF: Cartridge::chr_access<1>(addr, v); break;  // CHR-ROM/RAM
			case 0x2000 ... 0x3EFF: ciRam[nt_mirror(addr)] = v; break;         // Nametables
			case 0x3F00 ... 0x3FFF: // Palettes
				if ((addr & 0x13) == 0x10) addr &= ~0x10;
				cgRam[addr & 0x1F] = v; break;
		}
	}
	template <bool write> u8 access(u16 index, u8 v) { // Access PPU through registers
		static u8 res, buffer; // VRAM read buffer
		static bool latch;  // Detect second reading
		if (write) { // Write into register
			res = v;
			switch (index) {
				case 0:  ctrl.r = v; tAddr.nt = ctrl.nt; break;    // PPUCTRL   ($2000)
				case 1:  mask.r = v; break;                        // PPUMASK   ($2001)
				case 3:  oamAddr = v; break;                       // OAMADDR   ($2003)
				case 4:  oamMem[oamAddr++] = v; break;             // OAMDATA   ($2004)
				case 5:                                            // PPUSCROLL ($2005)
					if (!latch) { fX = v & 7; tAddr.cX = v >> 3; } // First write
					else  { tAddr.fY = v & 7; tAddr.cY = v >> 3; } // Second write
					latch = !latch; break;
				case 6:                                             // PPUADDR   ($2006)
					if (!latch) { tAddr.h = v & 0x3F; }             // First write
					else        { tAddr.l = v; vAddr.r = tAddr.r; } // Second write
					latch = !latch; break;
				case 7: wr(vAddr.addr, v); vAddr.addr += ctrl.incr ? 32 : 1; // PPUDATA ($2007)
			}
		}
		else { // Read from register
			switch (index) { // PPUSTATUS ($2002)				
				case 2:  res = (res & 0x1F) | status.r; status.vBlank = 0; latch = 0; break;
				case 4:  res = oamMem[oamAddr]; break; // OAMDATA ($2004)
				case 7:                                // PPUDATA ($2007)
					if (vAddr.addr <= 0x3EFF) res = buffer, buffer = rd(vAddr.addr);
					else res = buffer = rd(vAddr.addr);
					vAddr.addr += ctrl.incr ? 32 : 1;
			}
		}
		return res;
	}
	template u8 access<0>(u16, u8); template u8 access<1>(u16, u8);
	// Calculate graphics addresses 
	inline u16 nt_addr() { return 0x2000 | (vAddr.r & 0xFFF); }
	inline u16 at_addr() { return 0x23C0 | (vAddr.nt << 10) | ((vAddr.cY / 4) << 3) | (vAddr.cX / 4); }
	inline u16 bg_addr() { return (ctrl.bgTbl * 0x1000) + (nt * 16) + vAddr.fY; }
	// Increment the scroll by one pixel
	inline void h_scroll() { if (!rendering()) return; if (vAddr.cX == 31) vAddr.r ^= 0x41F; else vAddr.cX++; }
	inline void v_scroll() {
		if (!rendering()) return;
		if (vAddr.fY < 7) vAddr.fY++;
		else {
			vAddr.fY = 0;
			if (vAddr.cY == 31) vAddr.cY = 0;
			else if (vAddr.cY == 29) { vAddr.cY = 0; vAddr.nt ^= 0b10; }
			else vAddr.cY++;
		}
	}
	// Copy scrolling data from loopy T to loopy V
	inline void h_update() { if (!rendering()) return; vAddr.r = (vAddr.r & ~0x041F) | (tAddr.r & 0x041F); }
	inline void v_update() { if (!rendering()) return; vAddr.r = (vAddr.r & ~0x7BE0) | (tAddr.r & 0x7BE0); }
	inline void reload_shift() { // Put new data into the shift registers
		bgShiftL = (bgShiftL & 0xFF00) | bgL, bgShiftH = (bgShiftH & 0xFF00) | bgH;
		atLatchL = (at & 1), atLatchH = (at & 2);
	}
	void clear_oam() { // Clear secondary OAM
		for (int i = 0; i < 8; i++) {
			secOam[i].id = 64;
			secOam[i].y = secOam[i].tile = secOam[i].attr = secOam[i].x = 0xFF;
			secOam[i].dataL = secOam[i].dataH = 0;
		}
	}
	void eval_sprites() { // Fill secondary OAM with the sprite infos for the next scanline
		int n = 0;
		for (int i = 0; i < 64; i++) {
			int line = (scanline == 261 ? -1 : scanline) - oamMem[i*4 + 0];
			// If the sprite is in the scanline, copy its properties into secondary OAM
			if (line >= 0 and line < spr_height()) {
				secOam[n].id   = i;
				secOam[n].y    = oamMem[i*4 + 0];
				secOam[n].tile = oamMem[i*4 + 1];
				secOam[n].attr = oamMem[i*4 + 2];
				secOam[n].x    = oamMem[i*4 + 3];
				if (++n >= 8) { status.sprOvf = true; break; }
			}
		}
	}
	void load_sprites() { // Load the sprite info into primary OAM and fetch their tile data
		u16 addr;
		for (int i = 0; i < 8; i++)	{
			oam[i] = secOam[i]; // Copy secondary OAM into primary
			// Different address modes depending on the sprite height
			if (spr_height() == 16) addr = ((oam[i].tile & 1) * 0x1000) + ((oam[i].tile & ~1) * 16);
			else addr = (ctrl.sprTbl * 0x1000) + (oam[i].tile * 16);
			unsigned sprY = (scanline - oam[i].y) % spr_height();  // Line inside the sprite
			if (oam[i].attr & 0x80) sprY ^= spr_height() - 1;      // Vertical flip
			addr += sprY + (sprY & 8);  // Select the second tile if on 8x16
			oam[i].dataL = rd(addr + 0), oam[i].dataH = rd(addr + 8);
		}
	}
	void pixel() { // Process a pixel, draw it if it's on screen
		u8 palette = 0, objPalette = 0;
		bool objPriority = 0;
		int x = dot - 2;
		if (scanline < 240 and x >= 0 and x < 256) {
			if (mask.bg and not (!mask.bgLeft && x < 8)) { // Background
				palette = (NTH_BIT(bgShiftH, 15 - fX) << 1) | NTH_BIT(bgShiftL, 15 - fX);
				if (palette) palette |= ((NTH_BIT(atShiftH,  7 - fX) << 1) | NTH_BIT(atShiftL,  7 - fX))      << 2;
			}
			if (mask.spr and not (!mask.sprLeft && x < 8)) { // Sprites
				for (int i = 7; i >= 0; i--) {
					if (oam[i].id == 64) continue; // Void entry
					unsigned sprX = x - oam[i].x;
					if (sprX >= 8) continue; // Not in range
					if (oam[i].attr & 0x40) sprX ^= 7; // Horizontal flip
					u8 sprPalette = (NTH_BIT(oam[i].dataH, 7 - sprX) << 1) | NTH_BIT(oam[i].dataL, 7 - sprX);
					if (sprPalette == 0) continue; // Transparent pixel
					if (oam[i].id == 0 && palette && x != 255) status.sprHit = true;
					sprPalette |= (oam[i].attr & 3) << 2;
					objPalette = sprPalette + 16;
					objPriority = oam[i].attr & 0x20;
				}
			}
			if (objPalette && (palette == 0 || objPriority == 0)) palette = objPalette; // Evaluate priority
			pixels[scanline*256 + x] = nesRgb[rd(0x3F00 + (rendering() ? palette : 0))];
		}
		// Perform background shifts
		bgShiftL <<= 1, bgShiftH <<= 1;
		atShiftL = (atShiftL << 1) | atLatchL, atShiftH = (atShiftH << 1) | atLatchH;
	}
	template<Scanline s> void scanline_cycle() { // Execute a cycle of a scanline
		static u16 addr;
		if (s == NMI and dot == 1) { status.vBlank = true; if (ctrl.nmi) CPU::set_nmi(); }
		else if (s == POST and dot == 0) GUI::new_frame(pixels);
		else if (s == VISIBLE or s == PRE) {
			switch (dot) { // Sprites
				case   1: clear_oam(); if (s == PRE) { status.sprOvf = status.sprHit = false; } break;
				case 257: eval_sprites(); break;
				case 321: load_sprites(); break;
			}
			switch (dot) { // Background
				case 2 ... 255: case 322 ... 337:
					pixel();
					switch (dot % 8) {
						// Nametable
						case 1:  addr  = nt_addr(); reload_shift(); break;
						case 2:  nt    = rd(addr);  break;
						// Attribute
						case 3:  addr  = at_addr(); break;
						case 4:  at    = rd(addr);  if (vAddr.cY & 2) at >>= 4;
													if (vAddr.cX & 2) at >>= 2; break;
						// Background (low bits)
						case 5:  addr  = bg_addr(); break;
						case 6:  bgL   = rd(addr);  break;
						// Background (high bits)
						case 7:  addr += 8;         break;
						case 0:  bgH   = rd(addr); h_scroll(); break;
					} break;
				case         256:  pixel(); bgH = rd(addr); v_scroll(); break;  // Vertical bump
				case         257:  pixel(); reload_shift(); h_update(); break;  // Update horizontal position
				case 280 ... 304:  if (s == PRE)            v_update(); break;  // Update vertical position
				// No shift reloading
				case             1:  addr = nt_addr(); if (s == PRE) status.vBlank = false; break;
				case 321: case 339:  addr = nt_addr(); break;
				// Nametable fetch instead of attribute
				case           338:  nt = rd(addr); break;
				case           340:  nt = rd(addr); if (s == PRE && rendering() && frameOdd) dot++;
			}
			if (dot == 260 && rendering()) Cartridge::signal_scanline(); // Signal scanline to mapper
		}
	}
	void step() { // Execute a PPU cycle
		switch (scanline) {
			case 0 ... 239:  scanline_cycle<VISIBLE>(); break;
			case       240:  scanline_cycle<POST>();    break;
			case       241:  scanline_cycle<NMI>();     break;
			case       261:  scanline_cycle<PRE>();     break;
		}		
		if (++dot > 340) { dot %= 341; if (++scanline > 261) scanline = 0, frameOdd ^= 1; } // Update dot and scanline counters
	}
	void reset() { // Reset PPU
		frameOdd = false;
		scanline = dot = 0;
		ctrl.r = mask.r = status.r = 0;
		memset(pixels, 0x00, sizeof(pixels)), memset(ciRam,  0xFF, sizeof(ciRam)), memset(oamMem, 0x00, sizeof(oamMem));
	}
	void save() { // Save state
		memcpy(ss.ciRam, ciRam, sizeof ciRam), memcpy(ss.cgRam, cgRam, sizeof cgRam), memcpy(ss.oamMem, oamMem, sizeof oamMem);
		memcpy(ss.oam, oam, sizeof oam), memcpy(ss.secOam, secOam, sizeof secOam), memcpy(ss.pixels, pixels, sizeof pixels);
		ss.ctrl = ctrl, ss.mask = mask, ss.status = status;
	}
	void load() { // Load state
		memcpy(ciRam, ss.ciRam, sizeof ciRam), memcpy(cgRam, ss.cgRam, sizeof cgRam), memcpy(oamMem, ss.oamMem, sizeof oamMem);
		memcpy(oam, ss.oam, sizeof oam), memcpy(secOam, ss.secOam, sizeof secOam), memcpy(pixels, ss.pixels, sizeof pixels);
		ctrl = ss.ctrl, mask = ss.mask, status = ss.status;
	}
}
namespace Cartridge {
	Mapper* mapper = nullptr; // Mapper chip
	// PRG-ROM access
	template <bool wr> u8 access(u16 addr, u8 v) { return !wr ? mapper->read(addr) : mapper->write(addr, v); }
	template u8 access<0>(u16, u8); template u8 access<1>(u16, u8);
	// CHR-ROM/RAM access
	template <bool wr> u8 chr_access(u16 addr, u8 v) { return !wr ? mapper->chr_read(addr) : mapper->chr_write(addr, v); }
	template u8 chr_access<0>(u16, u8); template u8 chr_access<1>(u16, u8);
	void signal_scanline() { mapper->signal_scanline(); }
	void load(const char* fileName) { // Load the ROM from a file
		FILE* f = fopen(fileName, "rb");
		fseek(f, 0, SEEK_END);
		int size = ftell(f);
		fseek(f, 0, SEEK_SET);
		u8* rom = new u8[size];
		fread(rom, size, 1, f);
		fclose(f);
		int mapperNum = (rom[7] & 0xF0) | (rom[6] >> 4);
		if (mapper != nullptr) delete mapper;
		switch (mapperNum) {
			case 0: mapper = new Mapper0(rom); break;
			case 1: mapper = new Mapper1(rom); break;
			case 2: mapper = new Mapper2(rom); break;
			case 3: mapper = new Mapper3(rom); break;
			case 4: mapper = new Mapper4(rom); break;
			case 7: mapper = new Mapper7(rom); break;
		}
		CPU::power(), PPU::reset(), APU::reset();
	}
}
namespace Joypad {
	u8 joypad_bits[2]; // Joypad shift registers
	bool strobe; // Joypad strobe latch
	u8 read_state(int n) { // Read joypad state (NES register format)
		if (strobe) return 0x40 | (GUI::get_joypad_state(n) & 1); // When strobe is high, it keeps reading A
		// Get the status of a button and shift the register
		u8 j = 0x40 | (joypad_bits[n] & 1);
		joypad_bits[n] = 0x80 | (joypad_bits[n] >> 1);
		return j;
	}
	void write_strobe(bool v) {
		// Read the joypad data on strobe's transition 1 -> 0
		if (strobe and !v) for (int i = 0; i < 2; i++) joypad_bits[i] = GUI::get_joypad_state(i);
		strobe = v;
	}
}
namespace GUI {
	// SDL structures
	SDL_Window *window;
	SDL_Renderer *renderer;
	SDL_Texture *gameTexture, *saveGameTexture;
	u8 const *keys;
	Sound_Queue* soundQueue;
	SDL_Scancode KEY_A = SDL_SCANCODE_A, KEY_B = SDL_SCANCODE_S, KEY_SELECT = SDL_SCANCODE_SPACE, KEY_START = SDL_SCANCODE_RETURN;
	SDL_Scancode KEY_UP = SDL_SCANCODE_UP, KEY_DOWN = SDL_SCANCODE_DOWN, KEY_LEFT = SDL_SCANCODE_LEFT, KEY_RIGHT = SDL_SCANCODE_RIGHT;
	SDL_Scancode KEY_SAVE = SDL_SCANCODE_Q, KEY_LOAD = SDL_SCANCODE_W; // Saving and loading
	void init() { // Initialize GUI
		// Initialize graphics system
		SDL_Init(SDL_INIT_VIDEO | SDL_INIT_AUDIO | SDL_INIT_JOYSTICK);
		SDL_SetHint(SDL_HINT_RENDER_SCALE_QUALITY, "linear");
		APU::init();
		soundQueue = new Sound_Queue;
		soundQueue->init(96000);
		// Initialize graphics structures
		window = SDL_CreateWindow("BadNES", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, WIDTH, HEIGHT, 0);
		renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
		SDL_RenderSetLogicalSize(renderer, WIDTH, HEIGHT);
		gameTexture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_STREAMING, WIDTH, HEIGHT);
		keys = SDL_GetKeyboardState(0);		
	}
	u8 get_joypad_state(int n) { // Get the joypad state from SDL
		u8 j = 0;
		j |= keys[KEY_A] << 0; j |= keys[KEY_B] << 1; j |= keys[KEY_SELECT] << 2; j |= keys[KEY_START] << 3;
		j |= keys[KEY_UP] << 4; j |= keys[KEY_DOWN] << 5; j |= keys[KEY_LEFT] << 6; j |= keys[KEY_RIGHT] << 7;
		return j;
	}
	void new_frame(u32* pixels) { SDL_UpdateTexture(gameTexture, NULL, pixels, WIDTH * sizeof(u32)); } // Send the rendered frame to the GUI
	void new_samples(const blip_sample_t* samples, size_t count) { soundQueue->write(samples, count); }
	void render() { // Render the screen
		SDL_RenderClear(renderer);
		// Draw the NES screen
		SDL_RenderCopy(renderer, gameTexture, NULL, NULL);
		SDL_RenderPresent(renderer);
	}
	void save() { CPU::save(), PPU::save(), saveGameTexture = gameTexture; } // Save state
	void load() { CPU::load(), PPU::save(), gameTexture = saveGameTexture; } // Load state
	void run(const char* file) { // Run the emulator
		SDL_Event e;
		// Framerate control
		u32 frameStart, frameTime;
		const int FPS = 60, DELAY = 1000.0f / FPS;
		Cartridge::load(file);
		while (1) {			
			frameStart = SDL_GetTicks();
			// Handle events
			while (SDL_PollEvent(&e)) {
				if (e.type == SDL_QUIT) return;
				if (e.type == SDL_KEYDOWN) {
					if (keys[KEY_SAVE]) save();
					else if (keys[KEY_LOAD]) load();
				}
			}
			CPU::run_frame();
			render();
			// Wait to mantain framerate
			frameTime = SDL_GetTicks() - frameStart;
			if (frameTime < DELAY) SDL_Delay((int)(DELAY - frameTime));
		}
	}
}

int main(int argc, char *argv[]) { GUI::init(); GUI::run(argv[1]); }