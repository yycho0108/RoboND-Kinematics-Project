# Project: Kinematics Pick & Place

## Kinematic Analysis

![DH.png](figures/DH.png)

Accordingly, we find the DH parameters of the KUKA arm to be as follows:

&alpha;   | a      | d    | &theta;
:--------:|:------:|:----:|:-------:
0         | 0      | 0.75 | q&#8321;
-&pi;/2   | 0.35   | 0    | q&#8322;-&pi;/2
0         | 1.25   | 0    | q&#8323;
-&pi;/2   | -0.054 | 1.50 | q&#8324;
&pi;/2    | 0      | 0    | q&#8325;
-&pi;/2   | 0      | 0    | q&#8326;
0         | 0      | 0.303| 0

These values are defined [here](kuka_kin.py#112).

As a simple validation, I set the joint angles to zero and ran [fetch\_tf.py](./kuka_arm/scripts/fetch_tf.py) in order to see if the values from the tf transform matched.

Sample output:
```bash
# To reproduce the results:
# roslaunch kuka_arm forward_kinematics.launch
# rosrun kuka_arm fetch_tf.py

===
d1 : 0.75
a1 : 0.35
a2 : 1.25
d4 : 1.5
a3 : -0.054
d7 : 0.303
===
...
```

It is straightforward to see that the values exactly correspond to the parameters reported in the table.


## Project Implementation