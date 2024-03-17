# Enrico Avventi 2010

import numpy as np
import slycot

A = np.array([[    0 ,       1 ,  0 ,        0 ,  0 ,         0 ],
            [ -0.5 , -0.0002 ,  0 ,        0 ,  0 ,         0 ],
            [    0 ,       0 ,  0 ,        1 ,  0 ,         0 ],
            [    0 ,       0 , -1 , -0.00002 ,  0 ,         0 ],
            [    0 ,       0 ,  0 ,        0 ,  0 ,         1 ],
            [    0 ,       0 ,  0 ,        0 , -2 , -0.000002 ]])
B = np.array([[ 1 ],
            [ 0 ],
            [ 1 ],
            [ 0 ],
            [ 1 ],
            [ 0 ]])
C = np.array([[ 1 , 0 , 1 , 0 , 1 , 0 ]])
D = np.array([[ 0 ]])
out = slycot.ab13dd('C', 'I', 'N', 'D', 6, 1, 1, A, np.eye(6), B, C, D)
print('--- Example for ab13dd ---')
print('The L_infty norm of the system is')
print(out[0])
print('The peak frequency is')
print(out[1])