# Energy Based Stepdetection (EBS) 
This repository contains the implementation of the method described in the paper 'An Energy Based Scheme for Reconstruction of Piecewise Constant Signals underlying Molecular Machines'

[![Linux Build Status](https://travis-ci.org/qubit-ulm/ebs.svg?branch=master)](https://travis-ci.org/qubit-ulm/ebs)
[![Windows Build Status](https://ci.appveyor.com/api/projects/status/of2bssgx58e8qyi3?svg=true)](https://ci.appveyor.com/project/jrosskopf/ebs)

# Introduction 
Analyzing the physical and chemical properties of single DNA based molecular machines such as polymerases and helicases often necessitates to track motion on the length scale of base pairs. Although high resolution instruments have been developed that are capable to reach that limit, individual steps are often times hidden by a significant amount of noise which complicates data processing. 
The EBS project implements an effective algorithm which detects steps in a high bandwidth signal by minimizing energy functionals. In a first step an efficient convex denoising scheme is applied which allows compression to tupels of amplitudes and plateau lengths. Thus more sophisticated methods for assigning steps to the tupel data while accounting for prior information can be used. To this end we employed a combinatorial optimization algorithm formulated on a graph.

# Usage
After building EBS, the executables will reside in `build/bin/`. You can call them from there. The documentation below assumes, that the binaries are there.

## General remarks
EBS consists of four executables:

* `lambdaopt` using a heuristic to determining the optimal TVDN regularization parameter
* `denoising` to remove the noise by solving the TVDN problem
* `level_generator` to generate a file containing the set of levels to cluster on
* `graph_processing` to use an modified alpha-expansion to cluster the denoised levels

Each program has documentation which details how to use it, what each of the parameters are, and how to use them:

    $ ./bin/lambdaopt --help

If one wants to see more details of the inner workings of a program, adding the debug flag might increase the verbosity of the output:

    $ ./bin/lambdaopt --debug

The [Matrix Market File format](http://math.nist.gov/MatrixMarket/formats.html) as the default input output format. The format ASCII based, allows comment lines, which begin with a percent sign. We use the "array" format for general dense vectors. Details how to handle this format in python or matlab can be found in the particular demos.

If you find a bug in EBS, have a problem using it or have a question about the method in general, feel free to open a github issue. The development team will try to answer the problem or fix the issue in a timely fashion. 

## Determining the regularization parameter lambda
The first step in a typical processing chain of a piecewise constant singnal is to determine a reasonable choice for the regularization pparameter lambda. Typically this requires a lot of twiddling. The program `lambdaopt` implements the heuristic we proposed in the paper to chose this parameter automatically.

The usage is as simple as:

    $ ./bin/lambdaopt $NOISY_DATA

`$NOISY_DATA` is the path to the matrix market formated file containing the noisy input vector. The sole output of this command is a floating point value of the optimal lambda.

## Denoising the dataset
Having the optimal lambda at hand, the next step in the processing chain is solving the TVDN on the noisy input file:

    $ ./bin/denoising --lambda $LAMBDA_OPT $NOISY_DATA > $DENOISED_DATA

`$NOISY_DATA` is again the path to the noisy input vector. `$LAMBDA_OPT` contains the value for the regularization parameter. Typically this is the output of the `lambdaopt` program. The solution of the TVDN optimization problem is written to stdout by default, but a file destination can be chosen by supplying the `--output` parameter. In the example above the output to stdout is redirected to a file `$DENOISED_DATA`.

## Clustering to a set of predefined levels 
The prerequisit for clustering is having a noise free signal, which we assume in a matrix market file `$DENOISED_DATA`. The task for this step as described in the paper is to cluster or assign a level from a predefined set to each sample in the noise free vector. So as a second prerequisite we require a vector containing the level set. 
One option would be to generate a problem specific set of levels which incorporates prior knowledge about the steps to expect in the signal. Another option is to simply build a equidistant grid between the minimal and maximal value in `$DENOISED_DATA`. This is exactly what the `level_generator` program is good for:

    $ ./bin/level_generator --level-distance $DISTANCE $DENOISED_DATA > $LEVEL_DATA

The output, which is written to stdout by default, contains a vector with a grid with spacing of $DISTANCE. In the example the output is redirected to a file $LEVEL_DATA.

Now everything is at hand to start the clustering process:

    $ ./bin/graph_processing --input $DENOISED_DATA --levels $LEVEL_DATA > $CLUSTERED_DATA  

Where in the example the `$DENOISED_DATA` is the output of `denoising` and the `$LEVEL_DATA` is a vector containing the level set in matrix market format. The output is written to stdout by default and contains a vector of the same length as `$DENOISED_DATA`. 

The energy, which is minimized via graph cuts by the `graph_processing` program, consists of three components: A data term, a smoothing term and a step height prior term. The relative weight, between this terms, can be adjusted by setting the parameters `--rho-d` (default value 100), `--rho-s` (default value 10) and `--rho-p` (which is disabled by default). We found the default values to work well on our test datasets. But the values may need to be tweaked due to the specific problem.

To turn on the step height prior, the parameter `--rho-p` has to be chosen > 0. Futher the parameter `--prior-distance` has to be set to the distance of two adjacent steps, which should NOT be penalized.


# Demos
The source package contains demo code which demonstrates the usage of the above described programs form either [Matlab](http://www.mathworks.com/products/matlab/) or [Python](http://www.python.org). The Matlab demo is in the `matlab/demo.m` file. In this file the simulation code is used to create new test data for each run. The Python one in `python/demo.py`. Here we use pre-generated test data from the `noisy_data.mm` file in the same directory. The result of each demo run should be a plot which should look like the picture below:

![Output plot of the demos](https://github.com/qubit-ulm/ebs/blob/master/demo_plot_full.png)
![Takeout of the demo plot](https://github.com/qubit-ulm/ebs/blob/master/demo_plot_outtake.png)


# Building

## Requirements
EBS has the following requirements on Unix systems

    gcc >= 4.8 or
    Clang >= 3.3
    cmake >= 2.8.12 

on Windows systems 
    
    Visual Studio C++ >= 2012
    cmake >= 2.8.12

## Build Process
EBS uses CMake as a build system and allows several flexible build configuration options. One can consult any of numerous CMake tutorials for further documentation, but this tutorial should be enough to get EBS ready to use.

First clone the git repository using the command line git client or your favorite gui interface. By default this will create a directory `ebs` in the current working directory. You should change to this repository root, after the clone operation has finished.

    $ git clone --recursive https://github.com/qubit-ulm/ebs.git 
    $ cd ebs 

Then, make a build directory.  The directory can have any name, not just 'build', but 'build' is sufficient.

    $ mkdir build                                     
    $ cd build

The next step is to run CMake to configure the project.  Running CMake is the equivalent to running `./configure` with autotools. In this step CMake will also download the boost library, extract and bootstrap it. Boost is a dependency of EBS and is used for different purposes. All necessary parts of boost are statically linked, so a boost installation is not necessary on the computer.

    $ cmake ..

Once CMake has finished, the process building the executable depends on the opteration system you are on. On Unix systems the `make` program will start the build.

    $ make

On Windows systems you should run from the Visual Studio Command Prompt `msbuild`

    $ msbuild ebs.sln

This will build the denoising as well as the graph processing part of EBS.

If the build fails and you cannot figure out why, please file a ticket on the github page the EBS developers will quickly help you figure it out.


# Attribution

Please cite `An Energy Based Scheme for Reconstruction of Piecewise Constant Signals observed in the Movement of Molecular Machines, (2015)`
[http://arxiv.org/abs/1504.07873]. If you find this code useful in your research and add your project or publication to `the users list
<https://github.com/qubit-ulm/ebs/blob/master/docs/users.md>`_.

The BibTeX entry for our paper is::

    @article{ebs,
       author = {},
        title = {An Energy Based Scheme for Reconstruction of Piecewise Constant Signals observed in the Movement of Molecular Machines},
      journal = {Biophysical Journal},
         year = 2015,
       volume = ,
        pages = {},
       eprint = {},
          doi = {}
    }

# License
The EBS project is licensed to you under the Apache License, Version 2.0 (the "License"); you may not use the code except in compliance with the License. You may obtain a copy of the License at

[http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0)

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.
