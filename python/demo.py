import os
import sys
import subprocess
import scipy.io as sio
import matplotlib.pyplot as plt

def denoise_signal():
    cmd = ['lambdaopt', '--input', 'noisy_data.mm']
    out = subprocess.check_output(cmd)
    lambda_opt = float(out) 
   
    cmd = ['denoising', 
           '--lambda', str(lambda_opt), 
           '--input', 'noisy_data.mm', 
           '--output', 'denoised_data.mm'];
    subprocess.check_call(cmd);


def cluster_signal():
    distance = 0.2
    cmd = ['level_generator', 
           '--level-distance', str(distance), 
           '--input', 'denoised_data.mm', 
           '--output', 'level_data.mm'];
    subprocess.check_call(cmd)

    cmd = ['graph_processing', 
           '--levels', 'level_data.mm', 
           '--rho-d', str(0.01),
           '--rho-s',  str(0.1), 
           '--rho-p',  str(0.2), 
           '--prior-distance', str(distance), 
           '--input', 'denoised_data.mm', 
           '--output', 'clustered_data.mm'];
    subprocess.check_call(cmd)


def plot_result():
    time_data = sio.mmread('time_data.mm')

    plt.plot(time_data, sio.mmread('noisy_data.mm'), 
             color=[0.8,0.8,0.8] ,marker='.', markersize=1.0, linestyle='None',
             label='Noisy full bandwith data set')
    plt.plot(time_data, sio.mmread('pwcs_data.mm'), 
             color='r' ,marker='.', markersize=1.0, linestyle='None',
             label='Pure step signal')
    plt.plot(time_data, sio.mmread('denoised_data.mm'), 
             color='g' ,marker='.', markersize=1.0, linestyle='None',
             label='TVDN result')
    plt.plot(time_data, sio.mmread('clustered_data.mm'), 
             color='b' ,marker='.', markersize=1.0, linestyle='None',
             label='Clustered result')
    
    plt.title('Noisy, denoised and clustered data')
    plt.xlabel('time/s')
    plt.ylabel('bead motion/nm')
    plt.legend()

    plt.show()


if __name__ == "__main__":
    if os.name == 'nt':
        os.environ['PATH'] = ';'.join([os.getenv('PATH'), '..\\build\\bin'])
    else:
        os.environ['PATH'] = ':'.join([os.getenv('PATH'), '../build/bin'])

    sys.stdout.write('Denoising ... ')
    sys.stdout.flush()
    denoise_signal()
    print('done')

    sys.stdout.write('Clustering ... ')
    sys.stdout.flush()
    cluster_signal()
    print('done')

    print('Plotting ... ')
    plot_result()
    print('done')
