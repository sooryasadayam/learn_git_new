from SloppyCell.ReactionNetworks import *

expt = Experiment('expt1')
data = {'net1':{'data1': {0: (0, 0.037),
                          2: (0.16575, 0.025),
                          4: (0.43225, 0.033),
                          6: (0.48175, 0.035),
                          8: (0.46395, 0.0325),
                          10: (0.4081, 0.0255),
                          12: (0.37765, 0.0265),
                          14: (0.384, 0.0255),
                          16: (0.4208, 0.02),
                          18: (0.384, 0.02),
                          20: (0.4005, 0.024),
                          25: (0.3916, 0.026),
                          30: (0.4043, 0.027),
                          40: (0.2444, 0.0275),
                          50: (0.1391, 0.022),
                          60: (0.12765, 0.0355),
                          },
                'data2': {0: (1, 0.084),
                          2: (0.9275, 0.046),
                          4: (0.7923, 0.038),
                          6: (0.7778, 0.032),
                          8: (0.7053, 0.033),
                          10: (0.6522, 0.037),
                          12: (0.5894, 0.039),
                          14: (0.5894, 0.04),
                          16: (0.6377, 0.03),
                          18: (0.6425, 0.028),
                          20: (0.6908, 0.03),
                          25: (0.6908, 0.031),
                          30: (0.7585, 0.032),
                          40: (0.8068, 0.04),
                          50: (0.9275, 0.046),
                          60: (0.971, 0.082),
                          }
                }
        }
expt.set_data(data)
