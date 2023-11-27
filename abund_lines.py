abund_lines = {
    'C': {
        # 'CH': [4313-15-10+1, 4323+15+10],  # Use 30 AA window
        # 'C2': [4736-15-10, 5164-15-10],  # Use 30 AA window
        'CH': [4310],  # Use 20 AA window
        #'CH': [4280, 4313-5+1-5, 4323+5+4],  # Use 10 AA window
        'C2': [4736-5, 5120, 5164-5],  # Use 10 AA window
        # 'C13': [4215.0, 4251.5, ],
        'C13': [4190.0, 4201.20, 4207.8, 4217.7, 4221.8, 4230.2, 4237.6,
                4243.8, 4249.2, 4249.8, 4254.05]
    },
    'N': {
        # 'CN': [4170, 8200, 8325, 8375],  # Use 30 AA window
        # 'CN': [4170, 8325, 8375],  # Use 30 AA window
        # 'CN': [4170-10-20, 7920, 8040, 8325],  # Use 30 AA window
        # 'CN_blue': [4170-10-20],  # Use 30 AA window
        'CN_B': [4140.00, 4212.50],  # Use 30 AA window
        'CN_R': [7990, 8325],  # Use 30 AA window
        # 'CN': [4160, 4210-8, 8325, 8375],  # Use 30 AA window
    },
    'O': {
        'I': [6300.30, 6363.80, ],  # Forbidden lines
    },
    'Na': {
        'I': [5889.95, 5895.92, 6154.23, 6160.75, 8183.25, 8194.75],
    },
    'Mg': {
        # 'MgH': [5137., 5155., 5172.68, 5183.60, 5210.1]
        #'I': [3829.36, 3986.75, 4167.27, 4571.10, 4702.99,
        #      5172.68, 5183.60, 5528.40, 7691.00, 8736.02, 8806.76],
        'I': [3829.36, 3986.75, 4167.27, 4571.10, 4702.99,
              5172.68, 5183.60, 5528.40, 8806.76],
        # 'I': [4571.10, 4702.99,
        #       5172.68, 5183.60, 5528.40],
    },
    'Al': {
        'I': [3961.52, 3944.00, 6696.02, 6698.67, 7835., 7836., 8772., 8773.89],
        },
    'Si': {
        'I': [3905.52, 4102.94, 5690.43, 5793.07, 5948.54, 7405.77, 7409.08, 7415.95, 7423.50],
        },
    'K': {
        'I': [7664.91, 7698.97],
        },
    'Ca': {
        'I': [4226.73, 4283.01, 4318.65, 4425.44, 4434.96, 4435.69, 4454.78,
              4455.89, 5265.56, 5588.76, 5594.47, 6102.72, 6122.22, 6162.17,
              6169.06, 6439.07, 6455.60, 6471.66, 6493.78, 6499.65],
        'II': [8498.02, 8542.09, 8662.14],
    },
    'Sc': {
        'II': [4246.82, 4314.08, 4325.00, 4400.39, 4415.54, 5526.79, 5667.15, 5669.04, 6245.62,
               6279.74, 6309.90, 6604.60],
    },
    'Ti': {
        'I': [3989.76, 3998.64, 4008.93, 4533.25, 4534.78, 4535.57, 4656.47,
              4681.91, 4981.73, 4999.50, 5192.97, 5210.39, 6554.22, 6556.06,
              6599.11],
        'II': [3913.46, 4025.12, 4028.34, 4053.83, 4161.53, 4163.63, 4184.31,
               4290.22, 4300.05, 4330.72, 4337.91, 4394.06, 4395.03, 4395.84,
               4399.77, 4417.71, 4418.33, 4441.73, 4443.80, 4444.55, 4450.48,
               4464.45, 4468.52, 4470.85, 4501.27, 4529.48, 4533.96, 4563.77,
               4571.97, 4589.91, 4657.20, 4708.66, 4779.98, 4798.53, 4805.09,
               5129.16, 5188.69, 5226.54, 5336.79],
        },
    'V': {
        'I': [4379.23, 6216.37, 6274.66, 6285.16, 6292.82, 6357.29, 6531.42,
              6785.01],
        'II': [3951.96, 4005.70],
        },
    'Cr': {
        'I': [4254.33, 4274.80, 4289.72, 4337.57, 4616.14, 4626.19, 4646.15,
              4652.16, 5206.04, 5208.42, 5345.80, 5409.77],
        },
    'Mn': {
        'I': [4030.75, 4033.06, 4034.48, 4041.36, 4754.05, 4823.53],
        },
    'Fe': {
        # Look at Shetrone 1996 Deep mixing paper for more in the redish range
        'convol': [4383.55, 4859.74, 5012.07, 5110.41, 5166.28, 5171.60,
                   5194.94, 5232.94, 5269.54,  5371.49, 5397.13, 5405.77,
                   5429.70, 5434.52, 5446.92, 5497.52, 5506.78, 6191.56,
                   6494.98, 8688.62],
        'I': [3743.36, 3753.61, 3758.23, 3765.54, 3767.19, 3786.68, 3787.88,
              3805.34, 3815.84, 3825.88, 3827.82, 3840.44, 3841.05, 3845.17,
              3849.97, 3850.82, 3852.57, 3856.37, 3865.52, 3878.02, 3878.57,
              3886.28, 3887.05, 3895.66, 3899.71, 3902.95, 3917.18, 3920.26,
              3922.91, 3940.88, 3949.95, 3977.74, 4005.24, 4007.27, 4014.53,
              4021.87, 4032.63, 4045.81, 4058.22, 4063.59, 4067.27, 4071.74,
              4076.63, 4109.80, 4132.06, 4132.90, 4134.68, 4137.00, 4143.41,
              4143.87, 4147.67, 4152.17, 4153.90, 4154.50, 4154.81, 4156.80,
              4157.78, 4174.91, 4175.64, 4181.76, 4184.89, 4187.04, 4187.80,
              4191.43, 4195.33, 4196.21, 4199.10, 4202.03, 4216.18, 4222.21,
              4227.43, 4233.60, 4238.81, 4247.43, 4250.12, 4250.79, 4260.47,
              4271.15, 4271.76, 4282.40, 4325.76, 4337.05, 4352.73, 4375.93,
              4383.55, 4404.75, 4407.71, 4415.12, 4422.57, 4427.31, 4430.61,
              4442.34, 4443.19, 4447.72, 4454.38, 4459.12, 4461.65, 4466.55,
              4476.02, 4489.74, 4494.56, 4528.61, 4531.15, 4592.65, 4602.94,
              4632.91, 4647.43, 4691.41, 4707.27, 4733.59, 4736.77, 4859.74,
              4871.32, 4872.14, 4890.76, 4891.49, 4903.31, 4918.99, 4920.50,
              4938.81, 4939.69, 4994.13, 5001.87, 5012.07, 5041.07, 5041.76,
              5049.82, 5051.63, 5068.77, 5079.74, 5083.34, 5110.41, 5123.72,
              5127.36, 5142.93, 5150.84, 5151.91, 5166.28, 5171.60, 5191.45,
              5192.34, 5194.94, 5216.27, 5225.53, 5232.94, 5250.21, 5250.65,
              5254.96, 5266.56, 5269.54, 5281.79, 5283.62, 5302.30, 5307.36,
              5324.18, 5328.04, 5328.53, 5339.93, 5369.96, 5371.49, 5383.37,
              5393.17, 5397.13, 5405.77, 5415.20, 5429.70, 5434.52, 5446.92,
              5455.61, 5497.52, 5501.47, 5506.78, 5569.62, 5572.84, 5586.76,
              5615.64, 5624.54, 5658.82, 6065.48, 6136.61, 6137.69, 6191.56,
              6230.72, 6252.56, 6265.13, 6335.33, 6393.60, 6400.00, 6421.35,
              6430.85, 6494.98, 6592.91, 6677.99, 8514.07, 8582.26, 8688.62,
              8824.22],
        'II': [4178.86, 4233.17, 4416.82, 4491.41, 4508.28, 4515.34, 4520.22,
               4522.63, 4555.89, 4583.84, 4923.93, 5018.45, 5197.58, 5276.00],
        },
    'Co': {
        'I': [3845.47, 3873.12, 3881.87, 3995.31, 4020.90, 4110.53, 4118.77,
              4121.32],
    },
    'Ni': {
        'I': [3783.52, 3807.14, 3858.30, 5476.90, 6327.60, 6532.89, 6586.33,
              6643.60, 6767.80, 6772.36],
    },
    'Cu': {
        'I': [5105.53, 5782.13],
    },
    'Zn': {
        'I': [4722.15, 4810.53],
    },
    'Sr': {
        'II': [4077.71, 4215.52, 10036.70],
    },
    'Y': {
        'I': [6222.58, 6435.00],
        'II': [3774.33, 3788.69, 3818.34, 3950.35, 4398.01, 4682.33, 4854.87,
               4883.68, 4900.11],
    },
    'Zr': {
        'I': [6127.46, 6134.57, 6140.46, 6143.18],
        'II': [3836.76, 3998.97, 4149.20, 4161.21, 4208.98, 4496.98, 6114.85],
    },
    'Mo': {
        'I': [5570.44],
    },
    'Tc': {
        'I': [4262.27],
    },
    'Rb': {
        'I': [7800.28],
    },
    'Ba': {
        'II': [4554.03, 4934.09, 5853.69, 6141.73, 6496.90],
    },
    'La': {
        'II': [4077.00, 4086.71, 4123.23, 4429.88, 6390.48],
    },
    'Ce': {
        'II': [4083.22, 4137.64, 4222.60, 4418.78, 4562.36, 4628.16],
    },
    'Pr': {
        'II': [4062.80, 4408.84],
    },
    'Nd': {
        'II': [3784.25, 3838.98, 3991.74, 4012.70, 4021.33, 4023.00, 4061.08,
               4069.26, 4109.07, 4109.45, 4135.32, 4232.37, 4446.38, 4462.98,
               5319.81],
    },
    'Sm': {
        'II': [4318.93, 4433.80, 4434.32, 4467.34, 4519.63, 4537.94, 4577.69],
    },
    'Eu': {
        'II': [4129.70, 4205.00, 4435.60, 4522.50, 6645.00],
    },
    'Gd': {
        'II': [3796.38, 3916.51, 4191.07],
    },
    'Dy': {
        'II': [3869.86, 3944.68, 3996.69, 4077.90, 4103.31, 5169.69],
    },
}