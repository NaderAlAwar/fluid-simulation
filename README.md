# fluid-simulation

## Prerequisites

Install the following packages

```
sudo apt install libglfw3-dev
sudo apt install glew-utils libglew-dev libglm-dev
sudo apt install libxrandr-dev xorg-dev
```

## Building

```
cd src/
mkdir release
cd release
cmake .. -DCMAKE_BUILD_TYPE=RELEASE
```

## Running

To run fluid_spiral_*

```
./src/simulation -fluid ./simulations/fluid_spiral_xy.txt -timestep 0.01 -fluid_type water
```

To run fluid_drop.txt

```
./src/simulation -fluid ./simulations/fluid_drop.txt -timestep 0.001 -fluid_type water -enable_jets
```

To run fluid_dam.txt or fluid_stairs_3.txt

```
./src/simulation -fluid ./simulations/fluid_dam.txt -timestep 0.001 -fluid_type water
```

To run fluid_smoke.txt

```
./src/simulation -fluid ./simulations/fluid_drop.txt -timestep 0.001 -fluid_type smoke
```

To run fluid_lava.txt

```
./src/simulation -fluid ./simulations/fluid_lava.txt -timestep 0.001 -fluid_type lava
```

## Commands

At run-time, the simulation can be controlled with the following commands:

a: start/pause the animation

s: show/hide the surface (shows fluid instead of only particles)

v: show/hide velocity arrows

m: show/hide particles

p: show/hide pressure

d: toggle between showing different velocity arrows
