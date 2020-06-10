# skeletal-sex-estimation
A GNU Octave function for skeletal sex estimation

The function estimate_sex comprises a self-intuitive graphical user interface for estimating sex from the femur, tibia, and humerus. The available “long_bone_classifier.mat” contains all relevant classification parameters based on modern Greek reference skeletal collection known as the Athens Collection. A testing dataset “sampledata.csv” containing all relevant measurements from three long bones of one male and one female individual from the Athens Collection is also bundled with the function for demonstrating purposes.

Use the long-bone-diaphyseal-CSG-Toolkit to extract the necessary measurements for sex estimation based on the aforementioned bones. You may also utilize the inspect_CSG function available in the latest version of the CSG_Toolkit to aggregate all extracted skeletal properties in a csv file compliant with the format required by the estimate_sex function.
