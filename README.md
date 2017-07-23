# Izar

This small plugins is a collection of useful filters for post-processing CFD computations

- Cone slicing
- Sector duplication
- Computing aerodynamic variables from conservative values
- Changing frame (absolute/relative formulation)
- Generate vector from scalars

vtkSMPTools is used for multithreading. MPI is not supported right now.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

### Prerequisites

This software depends on :
 - cmake
 - Paraview (tested with 5.4 release). A local compilation is necessary to be able to compile this plugin
   [PARAVIEW_README](https://gitlab.kitware.com/paraview/paraview/blob/master/README.md)
 - YAML-cpp (tested with 0.5.3 version). This library needs boost library [YAML-cpp](https://github.com/jbeder/yaml-cpp)


### Installing

To install on Linux:

```
git clone Izar
```
```
cd Izar
```
```
mkdir build 
```
```
cd build 
```
```
ccmake ../src
```

Do the configuration to define ParaView compile location

```
make

```

You have a libizar.so that can be imported in ParaView plugin manager.

## Authors

* **Etienne Tang**
* **Mickael Philit** - [Mick7](https://github.com/MicK7)


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details




