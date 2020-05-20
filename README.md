# ABP-Durotaxis-Model
A simple model of active brownian motion used to verify the influence of durotaxis on cellular motility. 

## Background
Motivational work: https://doi.org/10.1103/PhysRevE.96.010402 

Much of the parameters used in this project are well defined in the above paper.

## Usage
From the working directory open and run any of the three ABP_Durotaxis scripts.

ABP_Durotaxis_Hist.py takes 7 input arguments (with inherent defaults if none are specified): 

  L: Transition Region Width, default = 3
  
  Dr: Effective Rotational Diffusion Constant, default = 5
  
  ks: soft stiffness parameter, default = 1
  
  kh: hard stiffness parameter, default = 50
  
  vel: cellular velocity, default = 1
  
  nw: number of cells used in simulation, default = 1000
  
  bb: bounding size of system, default = 20
  
  To run, execute in command shell: python ABP_Durotaxis_Hist.py L Dr ks kh vel nw bb (with actual values in place of the input parameters)
  
  The output of this script is a plot of cellular trajectories alongside a histogram of the final horizontal location data
  
ABP_Durotaxis_MSD.py takes 6 input parameters (w/ inherent defaults):

  L: Transition Region Width, default = 3
  
  ks: soft stiffness parameter, default = 1
  
  kh: hard stiffness parameter, default = 50
  
  dt: discrete timeset, default = 0.01
  
  nw: number of cells used in simulation, default = 1000
  
  ns: number of time steps used in simulation, default = 10^6
  
  To run, execute in command shell: python ABP_Durotaxis_MSD.py L ks kh dt nw ns (with actual values in place of input parameters)
  
  The output of this script is a loglog plot of mean square displacement of cells vs time
  
ABP_Durotaxis_DI.py takes 8 input parameters (w/ inherent defaults):

  L: Transition Region Width, default = 3
  
  Dr: Effective Rotational Diffusion Constant, default = 5
  
  ks: soft stiffness parameter, default = 1
  
  kh: hard stiffness parameter, default = 50
  
  vel: cellular velocity, default = 1
  
  nw: number of cells used in simulation, default = 1000
  
  bb: bounding size of system, default = 20
  
  ns: number of time steps, default = 1000
  
  To run, execute in command shell: python ABP_Durotaxis_MSD.py L Dr ks kh vel nw bb ns (w/ actual values in place of input parameters)
  
  The output of this script is a plot of durotaxis index (described in the paper in the background) vs time
  
  *The defaults for the parameters above are good starting points for running any of the scripts. There are a few guidelines for choosing parameter values if one wishes to stray from the defaults. L should be on the order of (and less than) bb, ks and kh should remain within an order of magnitude or two of each other (ks < kh), Dr and vel should be on the same order so that the cells exhibit APB. Other than that have fun playing around with the parameters.*
  
## Support
For support or questions reach out via nstclair@ucmerced.edu

## Acknowledgement
Thanks to Professor Kinjal Dabiswas of the department of Physics at UC Merced for getting the project off the ground and for the excellent guidance.

## License 
[MIT]

Copyright (c) [2020] [Nicholas St. Clair]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
