```
  ____            _ _   _ ______  _____ 
 |  _ \          | | \ | |  ____|/ ____|
 | |_) | __ _  __| |  \| | |__  | (___  
 |  _ < / _` |/ _` | . ` |  __|  \___ \ 
 | |_) | (_| | (_| | |\  | |____ ____) |
 |____/ \__,_|\__,_|_| \_|______|_____/ 
				EMULATOR
```

A nearly complete NES emulator in only 700 lines of C++, based on [LaiNES](https://github.com/AndreaOrru/LaiNES).

```
cloc main.cpp
       1 text file.
       1 unique file.                              
       0 files ignored.

github.com/AlDanial/cloc v 1.74  T=0.03 s (30.8 files/s, 24595.5 lines/s)
-------------------------------------------------------------------------------
Language                     files          blank        comment           code
-------------------------------------------------------------------------------
C++                              1              5             85            708
-------------------------------------------------------------------------------
```

## Features
- Cycle accurate
- Extremely compact: only 700 lines in one source file with no dependencies

## Requirements
- C++11
- SDL2

## Building
```sh
git clone --recursive https://github.com/Ta180m/BadNES && cd BadNES
g++ main.cpp -o badnes -std=c++11 -lSDL2main -lSDL2
```

## Usage
```sh
./badnes [path to rom]
```

## Compatibility
BadNES implements the most common mappers, which should be enough for a good percentage of the games:
- NROM (Mapper 000)
- MMC1 / SxROM (Mapper 001)
- UxROM (Mapper 002)
- CNROM (Mapper 003)
- MMC3, MMC6 / TxROM (Mapper 004)

You can check the compatibility for each ROM in the following list:
http://tuxnes.sourceforge.net/nesmapper.txt

## Credits
Special thanks to [Andrea Orru](https://github.com/AndreaOrru) for creating [LaiNES](https://github.com/AndreaOrru/LaiNES), the emulator that this project derives much of its code from.