- Current version: v0.7.1
- Completely tested version: v.0.4.0

--- Naming convention ---

- Single variable has 4-letter name, array has 5-letter, function has 6-letter

- Function naming:
  + ob*: obtain/observe
  + ad*: add
  + rm*: remove
  + mk*: make/set
  + ini*: initialize
  + nw*: update/calculate new
  + ap*: apply


--- I/O interface ---

- simula: Extracts data from source files to build the system, simulates the world.

- convrt: Converts the snapshot of the system at any moment from the output format to standard CHARMM format. Currently applied only the a system of water molecules.

- latGen: Creates lattice

- datAnl: Analyze data


--- File ---

- System class: world, tlist

- Object-classes: partk, bound, angle, dihed, resid (basically molecule, name chosen for convention)

- Type-classes (correspond to the object classes): ptype, btype, atype, dtype, rtype

- Memory-classes: memor, dmemo

- I/O interface: simula, convrt, latGen, datAnl

- Source file: param.txt, topol.txt, struc.txt, coord.txt, memor.txt


--- Data structure ---

- System-classes are overall containers
  + World: Contains every modules required for the simulation.
  + Tlist: Contains every modules required for the preparation of the simulation

- Type-classes are built during system construction, containing neccessary attributes to differentiate instances of object-classare. They are, however, not used for simulation: each instance of object-class will get and store all the necessary information required for simulation.
  + rtype: map associating input information with the respective type-classes.
  + ptype: name, charge, mass, epsilon and rmin: LJ-potential parameters (basically Van der Waals).
  + btype (atype, dtype): relaxed size, spring constant and pointers to ptype of the particles which make up the 2 (3, 4) ends of the bond (angle, dihedral angle).
    btype:  left --- righ
    atype:  left --- midl --- righ
    dtype:  uppr --- left --- righ --- lowr  

- Object-classes are the primary subject of the simulation:
  + Resid: Pointers to to objects (particle, bonds, angle) that make up the molecule.
  + Partk: Position, velocity, acceleration, list of bonds connecting to it, name.
  + Bound (angle, dihed): pointers to 2 (3, 4) partk that make up its ends and btype (atype, dtype).

- Memory-classes are just fancy cyclic queue to support MZ formalism


--- Source file ---

- param.txt: contains neccessary information to build the ptype, btype, atype, dtype

- topol.txt: contains neccessary information to build rtype

- struc.txt: contains neccessary information to build the structure of the simulated world

- coord.txt: contains the value of the initial time and neccessary information to initializes the position and velocity of the particles (every line corresponds to one particle, with the first 3 numbers are the components of the position vector, and the last 3 the velocity vector)

- memor.txt: contains the memory kernel (for now, it is simply the constant factor)


--- Output file ---

- trajc.trj: the output trajectory of the simulation

- coord.crd: the coordinate file in CHARMM-readable format
