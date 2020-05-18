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

// Initial declarations
namespace APU {
	template <bool write> u8 access(int elapsed, u16 addr, u8 v = 0);
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
	void init(), run();
	void new_frame(u32* pixels);
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
// PRG mapping functions
template <int pageKBs> void Mapper::map_prg(int slot, int bank) {
	if (bank < 0) bank = (prgSize / (0x400*pageKBs)) + bank;
	for (int i = 0; i < (pageKBs/8); i++) prgMap[(pageKBs/8) * slot + i] = (pageKBs*0x400*bank + 0x2000*i) % prgSize;
}
template void Mapper::map_prg<32>(int, int); template void Mapper::map_prg<16>(int, int); template void Mapper::map_prg<8> (int, int);
// CHR mapping functions
template <int pageKBs> void Mapper::map_chr(int slot, int bank) {
	for (int i = 0; i < pageKBs; i++) chrMap[pageKBs*slot + i] = (pageKBs*0x400*bank + 0x400*i) % chrSize;
}
template void Mapper::map_chr<8>(int, int); template void Mapper::map_chr<4>(int, int); template void Mapper::map_chr<2>(int, int); template void Mapper::map_chr<1>(int, int);
class Mapper0 : public Mapper {
	public: Mapper0(u8* rom) : Mapper(rom) { map_prg<32>(0, 0); map_chr<8> (0, 0); }
};
class Mapper1 : public Mapper {
	int writeN; u8 tmpReg, regs[4];
	// Apply the registers state 
	void apply() {
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
	// Apply the registers state
	void apply() {
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
	// Apply the registers state
	void apply() {
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
	// Apply the registers state
	void apply() {
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
	// Apply the registers state
    void apply() {
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
	const int OUT_SIZE = 4096;
	template <bool write> u8 access(int elapsed, u16 addr, u8 v) { return v; }
	template u8 access<0>(int, u16, u8); template u8 access<1>(int, u16, u8);
}
namespace CPU {
	u8 A, X, Y, S, ram[0x800]; u16 PC; Flags P; bool nmi, irq; // CPU state
	// Remaining clocks to end frame
	const int TOTAL_CYCLES = 29781; int remainingCycles;
	inline int elapsed() { return TOTAL_CYCLES - remainingCycles; }
	// Cycle emulation
	#define T   tick()
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
	// Flow control (branches, jumps)
	template<Flag f, bool v> void br() { 
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
	// Execute a CPU instruction
	void exec() {
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
	// Turn on the CPU
	void power() {
		remainingCycles = 0;
		P.set(0x04);
		A = X = Y = S = 0x00;
		memset(ram, 0xFF, sizeof(ram));
		nmi = irq = false;
		INT<RESET>();
	}
	// Run the CPU for roughly a frame
	void run_frame() {
		remainingCycles += TOTAL_CYCLES;
		while (remainingCycles > 0) {
			if (nmi) INT<NMI>();
			else if (irq and !P[I]) INT<IRQ>();
			exec();
		}
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
	u32 pixels[256 * 240];    // Video buffer
	Addr vAddr, tAddr; // Loopy V, T
	u8 fX, oamAddr; // Fine X, OAM address
	Ctrl ctrl;     // PPUCTRL   ($2000) register
	Mask mask;     // PPUMASK   ($2001) register
	Status status; // PPUSTATUS ($2002) register
	u8 nt, at, bgL, bgH, atShiftL, atShiftH; u16 bgShiftL, bgShiftH; // Background latches, shift registers
	bool atLatchL, atLatchH;
	int scanline, dot; bool frameOdd; // Rendering counters
	inline bool rendering() { return mask.bg || mask.spr; }
	inline int spr_height() { return ctrl.sprSz ? 16 : 8; }
	// Get CIRAM address according to mirroring
	u16 nt_mirror(u16 addr) {
		switch (mirroring) {
			case VERTICAL:   return addr % 0x800;
			case HORIZONTAL: return ((addr / 2) & 0x400) + (addr % 0x400);
			case ONE_SCREEN_LO:
        	case ONE_SCREEN_HI: return ((addr & 0x3ff) + ((mirroring == ONE_SCREEN_HI) ? 0x400 : 0x0)) - 0x2000;
			default:         return addr - 0x2000;
		}
	}
	void set_mirroring(Mirroring mode) { mirroring = mode; }
	// Access PPU memory
	u8 rd(u16 addr) {
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
	// Access PPU through registers
	template <bool write> u8 access(u16 index, u8 v) {
		static u8 res, buffer; // VRAM read buffer
		static bool latch;  // Detect second reading
		// Write into register
		if (write) {
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
	// Put new data into the shift registers
	inline void reload_shift() {
		bgShiftL = (bgShiftL & 0xFF00) | bgL, bgShiftH = (bgShiftH & 0xFF00) | bgH;
		atLatchL = (at & 1), atLatchH = (at & 2);
	}
	// Clear secondary OAM
	void clear_oam() {
		for (int i = 0; i < 8; i++) {
			secOam[i].id = 64;
			secOam[i].y = secOam[i].tile = secOam[i].attr = secOam[i].x = 0xFF;
			secOam[i].dataL = secOam[i].dataH = 0;
		}
	}
	// Fill secondary OAM with the sprite infos for the next scanline
	void eval_sprites() {
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
	// Load the sprite info into primary OAM and fetch their tile data
	void load_sprites() {
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
	// Process a pixel, draw it if it's on screen
	void pixel() {
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
	// Execute a cycle of a scanline
	template<Scanline s> void scanline_cycle() {
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
	// Execute a PPU cycle
	void step() {
		switch (scanline) {
			case 0 ... 239:  scanline_cycle<VISIBLE>(); break;
			case       240:  scanline_cycle<POST>();    break;
			case       241:  scanline_cycle<NMI>();     break;
			case       261:  scanline_cycle<PRE>();     break;
		}		
		if (++dot > 340) { dot %= 341; if (++scanline > 261) scanline = 0, frameOdd ^= 1; } // Update dot and scanline counters
	}
	void reset() {
		frameOdd = false;
		scanline = dot = 0;
		ctrl.r = mask.r = status.r = 0;
		memset(pixels, 0x00, sizeof(pixels)); memset(ciRam,  0xFF, sizeof(ciRam)); memset(oamMem, 0x00, sizeof(oamMem));
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
	// Load the ROM from a file.
	void load(const char* fileName) {
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
			case 0:  mapper = new Mapper0(rom); break;
			case 1:  mapper = new Mapper1(rom); break;
			case 2:  mapper = new Mapper2(rom); break;
			case 3:  mapper = new Mapper3(rom); break;
			case 4:  mapper = new Mapper4(rom); break;
			case 7:  mapper = new Mapper7(rom); break;
		}
		CPU::power(), PPU::reset();
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
	SDL_Window* window;
	SDL_Renderer* renderer;
	SDL_Texture* gameTexture;
	u8 const* keys;
	SDL_Scancode KEY_A = SDL_SCANCODE_A, KEY_B = SDL_SCANCODE_S, KEY_SELECT = SDL_SCANCODE_SPACE, KEY_START = SDL_SCANCODE_RETURN;
	SDL_Scancode KEY_UP = SDL_SCANCODE_UP, KEY_DOWN = SDL_SCANCODE_DOWN, KEY_LEFT = SDL_SCANCODE_LEFT, KEY_RIGHT = SDL_SCANCODE_RIGHT;
	// Initialize GUI
	void init() {
		// Initialize graphics system
		SDL_Init(SDL_INIT_VIDEO | SDL_INIT_AUDIO | SDL_INIT_JOYSTICK);
		SDL_SetHint(SDL_HINT_RENDER_SCALE_QUALITY, "linear");
		// Initialize graphics structures
		window = SDL_CreateWindow("BadNES", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, WIDTH, HEIGHT, 0);
		renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
		SDL_RenderSetLogicalSize(renderer, WIDTH, HEIGHT);
		gameTexture = SDL_CreateTexture (renderer, SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_STREAMING, WIDTH, HEIGHT);
		keys = SDL_GetKeyboardState(0);		
	}
	// Get the joypad state from SDL
	u8 get_joypad_state(int n) {
		u8 j = 0;
		j |= keys[KEY_A] << 0; j |= keys[KEY_B] << 1; j |= keys[KEY_SELECT] << 2; j |= keys[KEY_START] << 3;
		j |= keys[KEY_UP] << 4; j |= keys[KEY_DOWN] << 5; j |= keys[KEY_LEFT] << 6; j |= keys[KEY_RIGHT] << 7;
		return j;
	}
	// Send the rendered frame to the GUI
	void new_frame(u32* pixels) { SDL_UpdateTexture(gameTexture, NULL, pixels, WIDTH * sizeof(u32)); }
	// Render the screen
	void render() {
		SDL_RenderClear(renderer);
		// Draw the NES screen
		SDL_RenderCopy(renderer, gameTexture, NULL, NULL);
		SDL_RenderPresent(renderer);
	}
	// Run the emulator
	void run(const char* file) {
		SDL_Event e;
		// Framerate control
		u32 frameStart, frameTime;
		const int FPS   = 60;
		const int DELAY = 1000.0f / FPS;
		Cartridge::load(file);
		while (1) {			
			frameStart = SDL_GetTicks();
			// Handle events
			while (SDL_PollEvent(&e)) {
				if (e.type == SDL_QUIT) return;
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