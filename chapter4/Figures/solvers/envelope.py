
import numpy as np
import scipy.signal
import matplotlib.pyplot as plt
import scipy.interpolate

# all adapted from https://stackoverflow.com/questions/34235530/how-to-get-high-and-low-envelope-of-a-signal

def envelope(x,sig, distance=1):
    # split signal into negative and positive parts
    u_x = np.where(sig > 0)[0]
    l_x = np.where(sig < 0)[0]
    u_y = sig.copy()
    u_y[l_x] = 0
    l_y = -sig.copy()
    l_y[u_x] = 0
    
    # find upper and lower peaks
    u_peaks, _ = scipy.signal.find_peaks(u_y, distance=distance)
    l_peaks, _ = scipy.signal.find_peaks(l_y, distance=distance)
    
    # use peaks and peak values to make envelope
    u_x = u_peaks
    u_y = sig[u_peaks]
    l_x = l_peaks
    l_y = sig[l_peaks]
    
    # add start and end of signal to allow proper indexing
    end = len(sig)
    u_x = np.concatenate((u_x, [0, end]))
    u_y = np.concatenate((u_y, [0, 0]))
    l_x = np.concatenate((l_x, [0, end]))
    l_y = np.concatenate((l_y, [0, 0]))
    
    # create envelope functions
    u = scipy.interpolate.interp1d(u_x, u_y)
    l = scipy.interpolate.interp1d(l_x, l_y)
    return u, l

def envelope2(x,s):
    q_u = np.zeros(s.shape)
    q_l = np.zeros(s.shape)
    
    #Prepend the first value of (s) to the interpolating values. This forces the model to use the same starting point for both the upper and lower envelope models.
    
    u_x = [0,]
    u_y = [s[0],]
    
    l_x = [0,]
    l_y = [s[0],]
    
    #Detect peaks and troughs and mark their location in u_x,u_y,l_x,l_y respectively.
    
    for k in range(1,len(s)-1):
        if (np.sign(s[k]-s[k-1])==1) and (np.sign(s[k]-s[k+1])==1):
            u_x.append(k)
            u_y.append(s[k])
    
        if (np.sign(s[k]-s[k-1])==-1) and ((np.sign(s[k]-s[k+1]))==-1):
            l_x.append(k)
            l_y.append(s[k])
    
    #Append the last value of (s) to the interpolating values. This forces the model to use the same ending point for both the upper and lower envelope models.
    
    u_x.append(len(s)-1)
    u_y.append(s[-1])
    
    l_x.append(len(s)-1)
    l_y.append(s[-1])
    
    #Fit suitable models to the data. Here I am using cubic splines, similarly to the MATLAB example given in the question.
    
    u_p = scipy.interpolate.interp1d(u_x,u_y, kind = 'cubic',bounds_error = False, fill_value=0.0)
    l_p = scipy.interpolate.interp1d(l_x,l_y,kind = 'cubic',bounds_error = False, fill_value=0.0)
    #  u_p = scipy.interpolate.interp1d(u_x,u_y)
    #  l_p = scipy.interpolate.interp1d(l_x,l_y)
    
    #Evaluate each model over the domain of (s)
    #  for k in range(0,len(s)):
        #  q_u[k] = u_p(k)
        #  q_l[k] = l_p(k)
    return u_p, l_p
    
    
def envelope3(x,s, dmin=1, dmax=1, split=False):
    """
    Input :
    s: 1d-array, data signal from which to extract high and low envelopes
    dmin, dmax: int, optional, size of chunks, use this if the size of the input signal is too big
    split: bool, optional, if True, split the signal in half along its mean, might help to generate the envelope in some cases
    Output :
    lmin,lmax : high/low envelope idx of input signal s
    """

    # locals min      
    lmin = (np.diff(np.sign(np.diff(s))) > 0).nonzero()[0] + 1 
    # locals max
    lmax = (np.diff(np.sign(np.diff(s))) < 0).nonzero()[0] + 1 
    
    if len(lmin)==0 or len(lmax)==0:
        x=np.arange(len(s))
        se = scipy.interpolate.interp1d(x, s)
        return se,se
    if split:
        # s_mid is zero if s centered around x-axis or more generally mean of signal
        s_mid = np.mean(s) 
        # pre-sorting of locals min based on relative position with respect to s_mid 
        lmin = lmin[s[lmin]<s_mid]
        # pre-sorting of local max based on relative position with respect to s_mid 
        lmax = lmax[s[lmax]>s_mid]

    # global min of dmin-chunks of locals min 
    lmin = lmin[[i+np.argmin(s[lmin[i:i+dmin]]) for i in range(0,len(lmin),dmin)]]
    # global max of dmax-chunks of locals max 
    lmax = lmax[[i+np.argmax(s[lmax[i:i+dmax]]) for i in range(0,len(lmax),dmax)]]
    u_x=lmax
    u_y=s[lmax]
    l_x=lmin
    l_y=s[lmin]
    
    # add start and end of signal to allow proper indexing
    end = len(s)
    u_x = np.concatenate((u_x, [0, end]))
    u_y = np.concatenate((u_y, [u_y[0], u_y[-1]]))
    l_x = np.concatenate((l_x, [0, end]))
    l_y = np.concatenate((l_y, [l_y[0], l_y[-1]]))
    # create envelope functions
    u = scipy.interpolate.interp1d(u_x, u_y)
    l = scipy.interpolate.interp1d(l_x, l_y)
    return u, l




def envelope_monotonous(x,s, uvar="+", lvar="-"):
    
    if uvar=="+":
        u_x=[0]
        u_y=[s[0]]
        i=0
        while i<len(s)-1:
            i+=1
            if s[i]>u_y[-1] or i==len(s)-1:
                u_x.append(i)
                u_y.append(s[i])
    elif uvar=="-":
        u_x=[0]
        u_y=[max(s)]
        i=0
        while i<len(s)-1:
            xm=np.argmax(s[i+1:])+i+1
            u_x.append(x[xm])
            u_y.append(s[xm])
            i=xm
    else:
        raise ValueError("uvar should be + or -")
        
        
    if lvar=="-":
        l_x=[0]
        l_y=[s[0]]
        i=0
        while i<len(s)-1:
            i+=1
            if s[i]<u_y[-1] or i==len(s)-1:
                l_x.append(i)
                l_y.append(s[i])

    elif lvar=="+":
        l_x=[0]
        l_y=[min(s)]
        i=0
        while i<len(s)-1:
            xm=np.argmin(s[i+1:])+i+1
            l_x.append(x[xm])
            l_y.append(s[xm])
            i=xm
    else:
        raise ValueError("lvar should be + or -")
        


    

    

        #  print(i)

    # create envelope functions
    u = scipy.interpolate.interp1d(u_x, u_y)
    l = scipy.interpolate.interp1d(l_x, l_y)
    return u,l
def test(fun_envelope):
    x = np.arange(1000)
    sig = np.sin(x*25)* np.sin(x/2)*np.sin(x)/(x+1)**0.5
    u, l = fun_envelope(x,sig)
    
    plt.figure(figsize=(25,5))
    plt.plot(x, u(x), marker="x")
    plt.plot(x, l(x), marker="x")
    plt.plot(x, sig)


if __name__=="__main__":

    test(envelope)
    test(envelope2)
    test(envelope3)
    test(lambda x,s:envelope_monotonous(x,s,uvar="-",lvar="+"))
    plt.show()
