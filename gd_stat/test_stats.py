import numpy as np
from gd_stat import run_gd, stat_data, random_angles
from multiprocessing import Pool
import os
import time

# Generate seed from current time in milliseconds

########################################
#
# Computation Parameters
#
########################################
# angle_fn = 'actual_stat_angles'
# stat_fn = 'stats_data'
# eps = np.power(10, np.linspace(-5,-15,21, endpoint=True))
# angle_num = 5000
# seed = 42
angle_prefix = 'actual_stat_angles'
stat_fn = 'stats_data'
eps = np.power(10, np.linspace(-5, -15, 21, endpoint=True))
angle_num = 5000


def _one_run(args):
    j, ep, stat_fn, angle_prefix, angle_num = args

    data_fn = f'./{stat_fn}_{j}.txt'
    angle_fn = f'./{angle_prefix}_{j}.txt'

    # Unique seed per task (time + pid + j), still reproducible if you log it
    seed = (int(time.time_ns() // 1_000_000) ^ (os.getpid() << 16) ^ j) & 0xFFFFFFFF

    random_angles(size=angle_num, fn=angle_fn, seed=seed)

    run_gd(ep, theta_file=angle_fn, output_file=data_fn)

    data_temp = stat_data(data_file=data_fn, output_format='array').flatten()
    return j, seed, data_temp



if __name__ == '__main__':
    tasks = [(j, float(ep), stat_fn, angle_prefix, angle_num) for j, ep in enumerate(eps)]
    data_all = [None] * len(tasks)
    seeds = [None] * len(tasks)

    nproc = min(len(tasks), max(1, (os.cpu_count() or 1) - 1))
    with Pool(processes=nproc) as pool:
        results = pool.map(_one_run, tasks)

    for j, seed, data_temp in results:
        data_all[j] = data_temp
        seeds[j] = seed

    data_all = np.asarray(data_all)

    np.savetxt('all_stats.txt', data_all, fmt='%.16f')
    np.savetxt('all_stats_eps.txt', eps, fmt='%.16e')
    np.savetxt('all_stats_seeds.txt', np.asarray(seeds, dtype=np.uint32), fmt='%u')