import numpy as np
import sys
import os
import subprocess

GD_PATH = "../build/gd"

def random_angles(size:int = 10, fn:str = 'output.txt', seed:int=42):
    '''generate random angles, save them into a txt file'''

    # Create RNG with fixed seed for reproducibility
    rng = np.random.default_rng(seed)
    
    # Generate uniform random numbers in [0, 1)
    angles = rng.random(size) * 2 * np.pi
    
    # Save to text file
    np.savetxt(fn, angles, fmt='%.16f', delimiter='\n')
    
    return angles

def run_gd(epsilon:float, theta_file: str, output_file: str):
    '''run the gridsynth on the angle file with the given accuracy, then output the gate statistics in output files'''
    subprocess.run(
        [GD_PATH, f'{epsilon:.10e}', theta_file, output_file],
        check=True
    )

def stat_data(data_file: str, output_format = 'dict'):
    '''
    Get statistical data from CSV file.
    Report mean and std of the t, h, s, w gate counts.
    '''
    # Load data, skip header row
    data = np.loadtxt(data_file, delimiter=',', skiprows=1)
    
    # Columns:
    # 0 = theta
    # 1 = t_count
    # 2 = h_count
    # 3 = s_count
    # 4 = w_count
    
    t = data[:, 1]
    h = data[:, 2]
    s = data[:, 3]
    w = data[:, 4]
    
    if output_format == 'dict':
        stats = {
            "t_mean": np.mean(t),
            "t_std": np.std(t, ddof=1),
            "h_mean": np.mean(h),
            "h_std": np.std(h, ddof=1),
            "s_mean": np.mean(s),
            "s_std": np.std(s, ddof=1),
            "w_mean": np.mean(w),
            "w_std": np.std(w, ddof=1),
        }
    
        return stats
    else:
        stats = [[np.mean(t), np.std(t, ddof=1)],
        [np.mean(h), np.std(h, ddof=1)],
        [np.mean(w), np.std(s, ddof=1)],
        [np.mean(w), np.std(w, ddof=1)]]
        stats = np.array(stats)
        return stats

if __name__ == '__main__':
    random_angles(10, './angles.txt', seed=42)
    run_gd(1e-12, theta_file='./angles.txt', output_file='./output.txt')
    stats = stat_data('./output.txt', output_format='array')
    print('states: ')
    print(stats)