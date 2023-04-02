import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator


def convertDataFromTextToH5():
    nome_file = "hires_sym0_results" # from http://stellarcollapse.org/lippunerroberts2015
    ye_ar,s_ar,tau_ar,A,alpha,B1,beta1,B2,beta2,B3,beta3 = np.loadtxt(nome_file,
                                                                      unpack=True,
                                                                      usecols=(0,1,2,5,6,7,8,9,10,11,12))
    print(ye_ar.shape, len(set(ye_ar)), set(ye_ar))
    print(s_ar.shape, len(set(s_ar)), set(s_ar))
    print(tau_ar.shape, len(set(tau_ar)), set(tau_ar))
    print(A.shape, len(set(A)))
    print(alpha.shape, len(set(alpha)))
    print(B1.shape, len(set(B1)))

    # idx = np.where((s_ar == s1) & (tau_ar == tau1) & (ye_ar == ye1))

    new_ye = np.array(sorted(set(ye_ar)))
    new_s = np.array(sorted(set(s_ar)))
    new_tau = np.array(sorted(set(tau_ar)))
    new_A = []
    new_alpha = []
    new_B1 = []
    new_beta1 = []
    new_B2 = []
    new_beta2 = []
    new_B3 = []
    new_beta3 = []

    ii = 0
    for i, ye in enumerate(new_ye):
        for j, s in enumerate(new_s):
            for k, tau in enumerate(new_tau):
                idx = np.where((s_ar == s) & (tau_ar == tau) & (ye_ar == ye))
                new_A.append(A[idx])
                new_alpha.append(alpha[idx])
                new_B1.append(B1[idx])
                new_beta1.append(beta1[idx])
                new_B2.append(B2[idx])
                new_beta2.append(beta2[idx])
                new_B3.append(B3[idx])
                new_beta3.append(beta3[idx])
                # _idx = k + len(new_s) * (j + len(new_ye) * i)
                # print(f"i={i} j={j} k={k} ii={ii} idx={_idx}")
                ii = ii+1
                if (len(A[idx])>1):
                    raise ValueError()
    new_A = np.array(A)
    new_alpha = np.array(alpha)
    new_B1 = np.array(B1)
    new_beta1 = np.array(new_beta1)
    new_B2 = np.array(new_B2)
    new_beta2 = np.array(new_beta2)
    new_B3 = np.array(new_B3)
    new_beta3 = np.array(new_beta3)

    dfile = h5py.File("heatingLippuner.h5",'w')
    dfile.create_dataset(name="A",shape=new_A.shape, data=new_A)
    dfile.create_dataset(name="alpha",shape=new_alpha.shape, data=new_alpha)
    dfile.create_dataset(name="B1",shape=new_B1.shape, data=new_B1)
    dfile.create_dataset(name="beta1",shape=new_beta1.shape, data=new_beta1)
    dfile.create_dataset(name="B2",shape=new_B2.shape, data=new_B2)
    dfile.create_dataset(name="beta2",shape=new_beta2.shape, data=new_beta2)
    dfile.create_dataset(name="B3",shape=new_B3.shape, data=new_B3)
    dfile.create_dataset(name="beta3",shape=new_beta3.shape, data=new_beta3)

    dfile.create_dataset(name="Ye",shape=new_ye.shape, data=new_ye)
    dfile.create_dataset(name="s",shape=new_s.shape, data=new_s)
    dfile.create_dataset(name="tau",shape=new_tau.shape, data=new_tau)

    dfile.close()

    return (new_ye, new_s, new_tau, new_A, new_alpha, new_B1, new_beta1, new_B2, new_beta2, new_B3, new_beta3)

''' interpolate heating rate '''

def findIndex(x, arr, N):
    if x <= arr[0]:
        return 0
    elif x >= arr[N-1]:
        return N-2

    i = N // 2
    a = 0
    b = N-1

    # until the m_size of b-a > 1 continue shrinking the array, approaching the 'x'
    while b-a > 1:
        i = (b+a) // 2
        if arr[i] > x:
            b = i
        else:
            a = i

    return a

def intepPoint(p1, p2, p3):
    f3x = p1[3] + (p2[3] - p1[3]) * ((p3[0] - p1[0]) / (p2[0] - p1[0]))
    return f3x

def trilinear_interpolation(x1, y1, z1, f1, x2, y2, z2, f2, x3, y3, z3):
    # Calculate interpolation factors
    d = (x3 - x1) / (x2 - x1)
    e = (y3 - y1) / (y2 - y1)
    f = (z3 - z1) / (z2 - z1)

    # Get the values of f at the eight corners of the cube
    f3 = f1 + d * (f2 - f1)
    f4 = f1 + d * (f2 - f1) + e * (f2 - f1 - (f3 - f1))
    f5 = f1 + d * (f2 - f1) + (1 - e) * (f3 - f1)
    f6 = f1 + d * (f2 - f1) + e * (f2 - f1 - (f3 - f1)) + (1 - f) * (f3 - f1)
    f7 = f1 + (1 - d) * (f2 - f1) + e * (f2 - f1 - (f3 - f1))
    f8 = f1 + (1 - d) * (f2 - f1) + e * (f2 - f1 - (f3 - f1)) + (1 - f) * (f3 - f1)

    # Interpolate along the fourth dimension
    f3 = (1 - d) * ((1 - e) * ((1 - f) * f1 + f * f2) + e * ((1 - f) * f3 + f * f4)) + d * ((1 - e) * ((1 - f) * f5 + f * f6) + e * ((1 - f) * f7 + f * f8))

    return f3

def interp3D(ye_arr, s_arr, tau_arr, arr,
             ye_val, s_val, tau_val):

    ii = 0
    V = np.zeros((len(ye_arr),len(s_arr),len(tau_arr)))
    for i, ye in enumerate(ye_arr):
        for j, s in enumerate(s_arr):
            for k, tau in enumerate(tau_arr):
                V[i,j,k] = arr[ii]
                ii = ii + 1

    if (ye_val < ye_arr[0]):
        ye_val = ye_arr[0]
    if (ye_val > ye_arr[-1]):
        ye_val = ye_arr[-1]

    if (s_val < s_arr[0]):
        s_val = s_arr[0]
    if (s_val > s_arr[-1]):
        s_val = s_arr[-1]

    if (tau_val < tau_arr[0]):
        tau_val = tau_arr[0]
    if (tau_val > tau_arr[-1]):
        tau_val = tau_arr[-1]

    fn = RegularGridInterpolator((ye_arr,s_arr,tau_arr), V)
    val1 = np.float64( fn(np.array([ye_val,s_val,tau_val])) )
    print("{:.4e}".format(val1))

    _i_ye = findIndex(ye_val, ye_arr, len(ye_arr))
    _i_s = findIndex(s_val, s_arr, len(s_arr))
    _i_tau = findIndex(tau_val, tau_arr, len(tau_arr))

    _ip1_ye = _i_ye + 1
    if (_i_ye == len(ye_arr)-1):
        _ip1_ye = _i_ye

    _ip1_s = _i_s + 1
    if (_i_s == len(s_arr)-1):
        _ip1_s = _i_s

    _ip1_tau = _i_tau + 1
    if (_i_tau == len(tau_arr)-1):
        _ip1_tau = _i_tau

    # _idx = k + len(new_s) * (j + len(new_ye) * i)
    _i_val = arr[ _i_tau + len(s_arr) * (_i_s + len(ye_arr) * _i_ye) ]
    _ip1_val = arr[ _ip1_tau + len(s_arr) * (_ip1_s + len(ye_arr) * _ip1_ye) ]

    val = trilinear_interpolation(ye_arr[_i_ye], s_arr[_i_s], tau_arr[_i_tau], _i_val,
                                  ye_arr[_ip1_ye], s_arr[_ip1_s], tau_arr[_ip1_tau], _ip1_val,
                                  ye_val, s_val, tau_val)
    val = np.float64(val)
    print("my={:.4e} scipy={:.4e}".format(val,val1))
    # if (val < min(_i_val,_ip1_val) or val > max(_i_val, _ip1_val)):
    #     raise ValueError(f"val={val} < min(_i_val,_ip1_val)={min(_i_val,_ip1_val)} "
    #                      f"or val > max(_i_val, _ip1_val)={max(_i_val, _ip1_val)}")
    return _ip1_val


def getHeatingRate(t_arr, ye_=0.21, s_=105.5, tau_=21.5):

    dfile = h5py.File("heatingLippuner.h5",'r')
    new_ye = np.array(dfile["Ye"])
    new_s = np.array(dfile["s"])
    new_tau = np.array(dfile["tau"])
    # ----------------------------------
    new_alpha = np.array(dfile["alpha"])
    new_A = np.array(dfile["A"])
    new_B1 = np.array(dfile["B1"])
    new_beta1 = np.array(dfile["beta1"])
    new_B2 = np.array(dfile["B2"])
    new_beta2 = np.array(dfile["beta2"])
    new_B3 = np.array(dfile["B3"])
    new_beta3 = np.array(dfile["beta3"])


    # Flat[x + WIDTH * (y + DEPTH * z)] = Original[x, y, z]
    A_ = interp3D(new_ye, new_s, new_tau, new_A, ye_, s_, tau_)
    alpha_ = interp3D(new_ye, new_s, new_tau, new_alpha, ye_, s_, tau_)
    B1_ = interp3D(new_ye, new_s, new_tau, new_B1, ye_, s_, tau_)
    beta1_ = interp3D(new_ye, new_s, new_tau, new_beta1, ye_, s_, tau_)
    B2_ = interp3D(new_ye, new_s, new_tau, new_B2, ye_, s_, tau_)
    beta2_ = interp3D(new_ye, new_s, new_tau, new_beta2, ye_, s_, tau_)
    B3_ = interp3D(new_ye, new_s, new_tau, new_B3, ye_, s_, tau_)
    beta3_ = interp3D(new_ye, new_s, new_tau, new_beta3, ye_, s_, tau_)


    Q_arr = []
    for i in range(len(t_arr)):
        Q = A_ * pow(t_arr[i], -alpha_) \
          + B1_ * np.exp(-t_arr[i] / beta1_) \
          + B2_ * np.exp(-t_arr[i] / beta2_) \
          + B3_ * np.exp(-t_arr[i] / beta3_)
        Q_arr.append(np.log10(Q))

    # print(A_, new_A.min(), new_A.max())
    return np.array(Q_arr)

def smoothclamp(x, x1, x2, y1, y2):
    return y1 + (y2-y1)*(lambda t: np.where(t < 0 , 0, np.where( t <= 1 , 3.*t**2-2.*t**3, 1 ) ) )( np.log10(x/x1)/np.log10(x2/x1) )

def yevar_eps_nuc(a_eps_nuc,b_eps_nuc,t_eps_nuc,time):
    time_day=time*1.157407407e-5 #[day/s]
    if (np.isscalar(time_day)):
        tmp = min(max(4*time_day-4.,-20),20)   # t_eps_nuc still missing!
    else:
        T1=np.zeros(len(time_day))
        T1[(4.*time_day-4.)>(-20)]=4.*time_day[(4.*time_day-4.)>(-20)]-4.
        T1[(4.*time_day-4.)<(-20)]=-20.
        tmp=np.zeros(len(time_day))
        tmp[T1<20.]=T1[T1<20.]
        tmp[T1>20.]=20.
    return a_eps_nuc + b_eps_nuc/(1.+np.exp(tmp))

def calc_eps_nuc(kappa,time,eps0,a_eps_nuc,b_eps_nuc,t_eps_nuc):
    weight=smoothclamp(kappa,1.,10.,1.,0.)
    return eps0 * ( (1.- weight) + weight * yevar_eps_nuc(a_eps_nuc,b_eps_nuc,t_eps_nuc,time) )


def main():
    ye_arr = np.asarray([0.4848837209302326,0.4965116279069767,0.4906976744186047])
    s_arr = np.asarray([1.0,1.3,1.8])
    v_arr = np.asarray([3.e+9,2.e+9,1.e+9])
    t_array=np.logspace(-2.0, 2.0, num=50, endpoint=True)
    fig, ax1 = plt.subplots(1,1)
    for ye, s, v in zip(ye_arr, s_arr, v_arr):
        radius=1.0E+8 #cm
        tau=radius / v*1.0e+3 #tau e' ms, v in cm s^-1
        heat_rate = getHeatingRate( t_array, ye, s, tau )
        ax1.loglog(t_array,10**heat_rate,'-')
    plt.tight_layout()
    # plt.savefig('heating.pdf')
    plt.show()
    #  Lippuner+ 2015
    # nome_file = "hires_sym0_results"
    # ye_ar,s_ar,tau_ar,A,alpha,B1,beta1,B2,beta2,B3,beta3 = np.loadtxt(nome_file,
    #                                                                   unpack=True,
    #                                                                   usecols=(0,1,2,5,6,7,8,9,10,11,12))
    # print(ye_ar.shape, len(set(ye_ar)), set(ye_ar))
    # print(s_ar.shape, len(set(s_ar)), set(s_ar))
    # print(tau_ar.shape, len(set(tau_ar)), set(tau_ar))
    # print(A.shape, len(set(A)))
    # print(alpha.shape, len(set(alpha)))
    # print(B1.shape, len(set(B1)))
    #
    # ye1 = 0.22
    # s1 = 100.0
    # tau1 = 21.0
    # idx = np.where((s_ar == s1) & (tau_ar == tau1) & (ye_ar == ye1))
    #
    # new_ye = np.array(sorted(set(ye_ar)))
    # new_s = np.array(sorted(set(s_ar)))
    # new_tau = np.array(sorted(set(tau_ar)))
    # new_A = []
    # new_alpha = []
    # new_B1 = []
    # new_beta1 = []
    # new_B2 = []
    # new_beta2 = []
    # new_B3 = []
    # new_beta3 = []
    #
    #
    # for i, ye in enumerate(new_ye):
    #     for j, s in enumerate(new_s):
    #         for k, tau in enumerate(new_tau):
    #             idx = np.where((s_ar == s) & (tau_ar == tau) & (ye_ar == ye))
    #             new_A.append(A[idx])
    #             new_alpha.append(alpha[idx])
    #             new_B1.append(B1[idx])
    #             new_beta1.append(new_beta1[idx])
    #             new_B2.append(new_B2[idx])
    #             new_beta2.append(new_beta2[idx])
    #             new_B3.append(new_B3[idx])
    #             new_beta3.append(new_beta3[idx])
    #
    #             if (len(A[idx])>1):
    #                 raise ValueError()
    #
    #
    # new_A = np.array(A[idx])
    # new_alpha = np.array(alpha[idx])
    # new_B1 = np.array(B1[idx])
    # new_beta1 = np.array(new_beta1[idx])
    # new_B2 = np.array(new_B2[idx])
    # new_beta2 = np.array(new_beta2[idx])
    # new_B3 = np.array(new_B3[idx])
    # new_beta3 = np.array(new_beta3[idx])
    #




    print("Done")

# # Define the two points
# p1 = [1, 2, 3, 4]
# p2 = [5, 6, 7, 8]
#
# # Define the point to interpolate
# p3 = [2, 4, 5]
#
# # Perform linear interpolation
# f3x = p1[3] + (p2[3] - p1[3]) * ((p3[0] - p1[0]) / (p2[0] - p1[0]))
# f3y = p1[3] + (p2[3] - p1[3]) * ((p3[1] - p1[1]) / (p2[1] - p1[1]))
# f3z = p1[3] + (p2[3] - p1[3]) * ((p3[2] - p1[2]) / (p2[2] - p1[2]))
#
# # Plot the results
# fig, axs = plt.subplots(3, 1, figsize=(6, 10))
#
# axs[0].scatter([p1[0], p2[0], p3[0]], [p1[3], p2[3], f3x])
# axs[0].set_xlabel('x')
# axs[0].set_ylabel('f')
#
# axs[1].scatter([p1[1], p2[1], p3[1]], [p1[3], p2[3], f3y])
# axs[1].set_xlabel('y')
# axs[1].set_ylabel('f')
#
# axs[2].scatter([p1[2], p2[2], p3[2]], [p1[3], p2[3], f3z])
# axs[2].set_xlabel('z')
# axs[2].set_ylabel('f')
#
# plt.tight_layout()
# plt.show()

if __name__ == '__main__':
    main()