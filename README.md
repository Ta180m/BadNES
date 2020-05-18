```
  ____            _ _   _ ______  _____ 
 |  _ \          | | \ | |  ____|/ ____|
 | |_) | __ _  __| |  \| | |__  | (___  
 |  _ < / _` |/ _` | . ` |  __|  \___ \ 
 | |_) | (_| | (_| | |\  | |____ ____) |
 |____/ \__,_|\__,_|_| \_|______|_____/ 
				EMULATOR
```

A nearly complete NES emulator in a single 700-line C++ source file, based on [LaiNES](https://github.com/AndreaOrru/LaiNES).

```
cloc main.cpp
       1 text file.
       1 unique file.                              
       0 files ignored.

github.com/AlDanial/cloc v 1.74  T=0.04 s (25.1 files/s, 20506.9 lines/s)
-------------------------------------------------------------------------------
Language                     files          blank        comment           code
-------------------------------------------------------------------------------
C++                              1              5             86            727
-------------------------------------------------------------------------------
```

## Features
- Cycle accurate
- Extremely compact: only 700 lines in one source file with no dependencies

## Requirements
- C++11
- SDL2

## Building
### *nix systems
```sh
git clone --recursive https://github.com/Ta180m/BadNES && cd BadNES
g++ main.cpp -o badnes -std=c++11 -lSDL2main -lSDL2 -O3
```

### Windows
```sh
git clone --recursive https://github.com/Ta180m/BadNES && cd BadNES
g++ main.cpp -o badnes -std=c++11 -IC:\mingw\include\SDL2 -LC:\mingw\lib -w -Wl,-subsystem,windows -lmingw32 -lSDL2main -lSDL2 -O3
```

## Usage
### *nix systems
```sh
./badnes [path to ROM]
```

### Windows
Drag ROM over `badnes.exe`

## Compatibility
BadNES implements the most common mappers, which should be enough for a good percentage of the games:
- NROM (Mapper 000)
- MMC1 / SxROM (Mapper 001)
- UxROM (Mapper 002)
- CNROM (Mapper 003)
- MMC3, MMC6 / TxROM (Mapper 004)
- AxROM (Mapper 007) (Not well supported, expect bugs!)

You can check the compatibility for each ROM in the following list:
http://tuxnes.sourceforge.net/nesmapper.txt

## FAQ
### Why the name?
http://www.usaco.org/index.php?page=viewproblem2&cpid=1041

### Why did you do this?
I really like LaiNES: it's very minimal and compact — but it has several dependencies, so I try to write the shortest NES emulator possible, in a single source file.

### There's no sound!
I'm working on it right now.

### What do you mean by "nearly complete"?
I'll say it's complete once I finish implementing the APU.

### This game doesn't work!
Check to make sure BadNES implements its mapper. Also, Mapper 7 is known to be especially buggy.

## Credits
Special thanks to [Andrea Orru](https://github.com/AndreaOrru) for creating [LaiNES](https://github.com/AndreaOrru/LaiNES), the emulator that this project derives much of its code from.