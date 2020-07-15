# skeletal-sex-estimation
A GNU Octave function for skeletal sex estimation

The function estimate_sex provides a graphical user interface for estimating sex from various sleletal elements according to the selected classifier saved in the respective .mat files. The required classifier file can be parsed as an input argument to the function or it can be selected through the open file dialog when the function is called withough any input arguments.
For further details regarding its usage type:

		help estimate_sex

Currently, the function supports two types of classifiers.

1) Use "long_bone_classifier.mat" for estimating sex from CSG data of the femur, tibia, and humerus bones. The CSG data can be automatically extracted from virtual 3D bone models with the long-bone-diaphyseal-CSG-Toolkit, which also saves the required data to the appropriate format.

2) Use "vertebrae_classifier.mat" for estimating sex based on linear measurements from the T1, T12 and L1 vertebrae.
