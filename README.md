# Unscented Kalman Filter
A Project of the Udacity Self-Driving Car Engineer Nanodegree Program
---

Kalman filter is a enginereing technique to fuse the different information to predict the next state of a subject of matter in changing environment.
Unscented kalman filter in particular improves the performance of the prediction by mitigating noises produces non-linearlity in the process.
In this proejct given the sample data from LIDAR/RADAR sensors Kalman filter run the cycle of preduction and state update as a simutation.
The evaluation of the performance is measured in Root of Mean Squared Error between the prediction and the actual.

## Dependencies

* cmake >= v3.5
* make >= v4.1
* gcc/g++ >= v5.4

### How to run and see result: 
After configurint the dependencies, please clone this repo and run `./build/UnscentedKF path/to/obj_pose-laser-radar-synthetic-input.txt path/to/output.txt`.