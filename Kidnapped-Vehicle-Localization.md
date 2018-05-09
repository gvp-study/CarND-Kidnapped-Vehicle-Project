
---

**Kidnapped Vehicle Localization Project**

The goal of this project are the following:

Implement the Particle Filter Localization method for tracking the pose of a car in a map with prediction of its state from its motion sensors (velocity, yaw_rate) and update of its state using measurements of landmarks from a local range sensor.

[//]: # (Image References)
[image1]: ./examples/first.png
[image2]: ./examples/second.png
[image3]: ./examples/third.png
[image4]: ./examples/fourth.png
[video1]: ,/examples/localization.mp4

## [Rubric Points](https://review.udacity.com/#!/rubrics/747/view)
All the code for this project has been derived from the example code in the course and is in this directory.
[Here](https://github.com/gvp-study/CarND-Kidnapped-Vehicle-Project.git)
---

### Kidnapped Vehicle Localization
## Implementing the Particle Filter
The directory structure of this repository is as follows:

```
root
|   build.sh
|   clean.sh
|   CMakeLists.txt
|   README.md
|   run.sh
|
|___data
|   |   
|   |   map_data.txt
|   
|   
|___src
|    |   helper_functions.h
|    |   main.cpp
|    |   map.h
|    |   particle_filter.cpp
|    |   particle_filter.h
|
|___examples
|    | images and videos   
```
I modified the four main functions (init, prediction updateWeights and resample) in particle_filter.cpp. These functions implement the main methods of the particle filter algorithm.
1. init: This function initializes the particle filter with a certain number of particles and sets them to cluster around the given mean pose with a standard deviations using a normal distribution. The probability of each particle at init is set to 1.0.
2. prediction: This function takes as input the elapsed time and the velocity and yaw_rate and predicts the pose of each particle. This prediction is computed using a normal noise distribution given by the standard deviation.
3. updateWeights: This function takes as input the set of observations from the range sensor and matches them to the nearest landmark for each particle. The resulting probability is stored as the updated weight of that particle.
4. resample: This function resamples the particles based on their probability (weight) and generates a new set of particles. This new set replaces the current samples.

After the initialization, the filter estimates the pose of the car first by prediction using the motion and then updating the pose based on the range sensor measurements and the landmarks in the map. The result is shown in the output of the simulator as shown in the figures below.

![alt text][image1]

![alt text][image2]

![alt text][image3]

The final pose and the error show that the implementaion was successful in terms of accuracy and speed.

![alt text][image4]

The movie of the simulator running the particle filter localization code is shown below.[link to my video](./examples/localization.mp4)
![alt text][video1]
The things the grading code is looking for are:


1. **Accuracy**: The max translation error and max yaw error are within bounds.

2. **Performance**: The filter completed within the time of 100 seconds.
