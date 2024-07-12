# Caustic-Design-Project
Based on the original paper: [High-Contrast Computational Caustic Design](https://taiya.github.io/pubs/schwartzburg2014caustics.pdf)

There are two programs in this project.

 1.  Caustic-Design - handles Optimal-Transport and corresponding tasks
 2.  Target Optimization - handles 3D Optimization

## Goal
The goal of this project is to create an open source application that enables artists and researchers to compute the geometry of a surface such that its caustics cast a specified target distribution (an image for example).

![animation](./animation.gif)

## Building

### Dependencies:

1. `sudo apt update`

2. Install Qt4:

```sh
sudo apt install libqt4-dev
```

If your system can not install `libqt4-dev` directly, you can try to install it from `ppa:rock-core/qt4` by adding the path first:
*(Tested on Debian10 and Ubuntu 20.04. Because use of `ppa:rock-core/qt4`, it is not compatible with Ubuntu 22.04 and Debian11.)*

```sh
sudo add-apt-repository ppa:rock-core/qt4
```

Then retry to install `libqt4-dev`:

```sh
audo apt install libqt4-dev
```

3. Install other dependencies:

```sh
sudo apt install cmake libblas-dev liblapack-dev libtbb-dev libmetis-dev build-essential libsuitesparse-dev liblbfgs-dev libtinyxml-dev libgmp3-dev libmpc-dev libboost-all-dev libglew-dev libsoil-dev libassimp-dev checkinstall
```

4. Install CGAL 4.9.1:

```sh
mkdir ~/source
cd ~/source
wget https://github.com/CGAL/cgal/releases/download/releases%2FCGAL-4.9.1/CGAL-4.9.1.tar.xz
tar -xf CGAL-4.9.1.tar.xz
mkdir build-CGAL-4.9.1
cd build-CGAL-4.9.1
cmake ../CGAL-4.9.1
make
sudo checkinstall
```

5. Install ceres-solver: http://ceres-solver.org/installation.html

```sh
cd ~/source
wget http://ceres-solver.org/ceres-solver-2.1.0.tar.gz
sudo apt update
sudo apt-get install libgoogle-glog-dev libgflags-dev libatlas-base-dev libeigen3-dev libsuitesparse-dev
tar zxf ceres-solver-2.1.0.tar.gz
mkdir ceres-bin
cd ceres-bin
cmake ../ceres-solver-2.1.0
make -j3
sudo make install
```

### Build
We suggest using cmake to build the project. To do so, simple:

Caustic Design:

```sh
mkdir build-Caustic_Design
cd build-Caustic_Design
cmake ../Caustic_Design/
make
```

Target Surface:

```sh
mkdir build-Target_Surface
cd build-Target_Surface
cmake ../Target_Surface/
make
```

Next, check for the existence of the dynamic library path environment variable (LD_LIBRARY_PATH)
`echo $LD_LIBRARY_PATH`

If there is nothing to be displayed, add a default path value
`export LD_LIBRARY_PATH=/usr/local/lib/`

If everything went well, you should be able to run Caustic_Design and Target_Surface by entering the command:
`./Caustic_Design`
or
`./Target_Surface`

## HOW TO USE
### Compute Optimal Transport
- First create your two grayscale images. The first image is the target image. This is the image that you want your lens or mirror to project. The second image is the source image. That is what the projection currently looks like which is in most cases a blank white image.
- Open Caustic_Design and click `File -> Load Image`, then load your target image. Next, click `File -> Load Source Image`, then load your source image.
- To set the optimal transport resolution, click `Algorithm -> Set Parameters` and change the `Sites` from 5000 to any other value. The higher this number is, the better your results will be, but the longer the optimal transport and the other steps will take. A site count of between 10k and 100k should be sufficient. You should also increase the `Levels` choose a value such that 4^(#levels) is about equal to the sites count. But ive noticed that the program crashes with a level count higher than 6, so i ket it at 6.
- Now finally click on `Algorithm -> Compute Optimal Transport`, click yes, and now the software will calculate the optimal transport between the source image and target image. It will take a couple of hours to complete depending on your specified steps and resolution. When it finishes computing, it will ask you to store the .weights file, save this file. You should also save the .dat file which you can do by clicking `File -> Save Source DAT`.

### Compute Interpolation
- Once you computed the optimal transport, the next thing to do is create a mesh in blender that will be deformed to create your caustic surface. I got the best results creating a cube, scaling it on the x-axis by 0.2, open it in edit mode, select the face facing the positive x-axis, click edge->subdivide, enter the amount of subdivisions you want, more results in higher resolution, but is slower to compute, then triangulate all faces by clicking 'Face' in the above menu and 'Triangulate Faces'.
- Then export this mesh as obj but make sure you unselect `Write Normals`, `Include UVs`, and `Write Materials` under `Geometry`.
- Next, launch the Target_Surface program and click `File -> Load Model`, then load your exported 3d model from blender, check to see if the subdivided face is facing the positive x-axis, and click `File -> Save Vertices`, save this file for the next step.
- Next, go back to Caustic_Design (this can be a new instance but doesnt have to be), Click on `Algorithm -> Compute Interpolation`, it will ask you for the target image again (the image you want to project), pointset (perviously generated .dat), weights (previously generated .weights), and light origin points (vertices generated from the previous step).
- when the last file is loaded, it will start the interpolation process. This will take a very long time depending on the resolution of the mesh and the amount of sites from the optimal transport step.
- When it is finished, it will ask you to save the interpolated point set (.dat), this file will be used to generate the caustic surface in the next step.

### Generate the Caustic Surface
- Open a new instance of the Target_Surface program and load your model again, then click on `File -> Load Light-Ray Receiver Position` open the .dat file from the interpolation step, if you then rotate the view, you should see your caustic image apear on the white target plane. if this is the case, click `Algorithm -> Run Target Optimization`. This will modify the surface of the model so it will project the caustic image. This will take a while depending on how easy it is to solve the surface hightmap.
- When the surface solver is done, you can export the model by clicking `File -> Save Model`, save the file as .stl
- You now have your caustic surface computed. This can be fabricated by a CNC machine. You can also simulate the caustics in blender using LuxRender.

## HOW DOES IT WORK
When we want to create a lens surface that casts a certain shadow, we want to change the surface as subtile as possible such that it redirects the light in just the right way to cast our disired image.

This is achieved in 3 steps:
- Get the current shadow of the lens
- find the optimal transport between the current shadow and the shadow we want to achieve
- encode this transport mapping in the surface of our lens

the last two steps are the most complex and fundamental to this method.

### Optimal transport
An optimal transport mapping is the most optimal way to move one distribution to another, in our case light.

the steps to compute the optimal transport are as follows:
- descretize the input distribution by dithering
- construct a power diagram (Laguerre–Voronoi diagram) of the descretization
- optimize the weights of the power diagram while continuing to update the power diagram so that it converges to a optimal mapping between the input and target distribution

In our case we dither the the input distribution using lloyd iterations on a power diagram with weights equal to 0 (regular voronoi diagram). The output of this descretization is then fed into the solver that finds the disired weights. 

when the solver finishes, we have a set of weights for the given dithered point set. When we construct the power diagram with the calculated weights, we get our target distribution, and if we construct the power diagram of the same points with weights equal to 0, we get our input distribution. This is our optimal transport mapping.

### Encoding the optimal transport on the surface
TODO

