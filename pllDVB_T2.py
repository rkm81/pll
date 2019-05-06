#!/usr/bin/env python
import time
import matplotlib.pyplot as plt
import random
import numpy as np

class registers_value(object):

	def __init__(self):
		self.soft_reload = 0

		self.gps_pll_lock = 0
		self.gps_holdover = 0
		self.pps_sync = 1

		self.gps_dac_value = 0x8000

		self.pd_error = 0
		self.pll_1pps_in_ok = 0
		self.pll_1pps_src = 0
		self.vco_10mhz_ctrl_src = 0
		self.pd_cntr_overflow = 0

		self.gps_user_dac_data = 0x8000

class bath(object):

	def __init__(self):
		self.fill_rate = 255
		self.drain_rate = 1
		self.thresh = 60e-9
		self.bathtub = 0 
		self.locked_level = 1024 
		self.unlocked_level = -1024
		self.lock  = 0

class filter(object):

	def __init__(self):
		self.y = [0, 0, 0, 0]
		self.x = [0, 0, 0]
		self.Kmax = 25
		self.Kmin = 8
		self.b = [-7.0344e-13 ,50.81332918, -49.97172275]
		self.a = [1963.636363, -5260.5274990, 4680.642293, -1383.751154]

def saturation(Xin, lower, upper):

	if Xin <= lower:
		Xin = lower
	elif Xin >= upper:
		Xin = upper
	return Xin

def filterFanction(Xin, gps_pll_lock , Filt):

	a0 = Filt.a[0]
	Filt.b = [x/a0 for x in Filt.b]	
	Filt.a = [s/a0 for s in Filt.a]	

	K = Filt.Kmax
	if gps_pll_lock  == 1:
		K = Filt.Kmin

	Filt.x.insert(0,Xin*K)
	Filt.x.pop()
	Filt.y.insert(0,0)
	Filt.y.pop()
	Filt.y[0] = Filt.b[0]*Filt.x[0] + Filt.b[1]*Filt.x[1] + Filt.b[2]*Filt.x[2] - Filt.a[1]*Filt.y[1] - Filt.a[2]*Filt.y[2] - Filt.a[3]*Filt.y[3]  
	Yout = Filt.y[0]
	return Yout

def PHD(t1, t2):

	Kphd = 100e6
	return saturation(round((t1 - t2) * Kphd), -2048, 2047)

def VCO(Yout, F0, Tout1):

	Kdac = 10/32768
	Kdiv = 1/10e6
	return ((saturation(round(Yout), -32768, 32767) * Kdac + F0) * Kdiv  + Tout1)

def lock_detector(bath_inf, err_in):

	if bath_inf.thresh >= abs(err_in):
		bucket = bath_inf.fill_rate
	else:
		bucket = -bath_inf.drain_rate
	
	bath_inf.bathtub = saturation(bath_inf.bathtub + bucket, -2048, 2048)

	if bath_inf.bathtub <= bath_inf.unlocked_level:
		bath_inf.lock  = 0
	elif bath_inf.bathtub >= bath_inf.locked_level:
		bath_inf.lock  = 1
	return bath_inf.lock 

def gist_PFD(Xin, G, j):

	numerator = 1
	mat = 0
	h = 0

	if j >=2:
		numerator = j - 1

	if Xin >= 101:
		G[len(G)-1] = G[len(G)-1] +1/numerator
	elif Xin <= -101:
		G[0] = G[0] +1/numerator
	else:
		G[101 + Xin] = G[101 + Xin] + 1/numerator

	for gg in range(0,203):
		G[gg] = G[gg] * (numerator/j)
		if gg >= 1 and gg <=201:
			mat = mat + (gg-100) * G[gg]
			h = h + ((gg-100)**2) * G[gg]

	sigm = round((h - (mat ** 2))**0.5)
	if sigm <= 2:
		sigm = 3

	return sigm

def sigma(buf_sigm, X):

	threshold = 100
	if len(buf_sigm) <= 1000:
		if X <= threshold and X >= -threshold:
			buf_sigm.append(X)

	else:
		if X <= threshold and X >= -threshold:
			buf_sigm.pop()
			buf_sigm.insert(0, X)

	# sigma = round(np.var(buf_sigm)**0.5)
	sigma = 20
	if sigma <= 2:
		sigma = 3

	return sigma

# def moving_average(X, buf):

# 	if len(buf) <= 500:
# 		buf.append(X)
# 	else:
# 		buf.pop()
# 		buf.insert(0,saturation(X, -32768, 32767))

# 	return round(sum(buf)/len(buf))

registers = registers_value()

t = []
out_err = []
input_err = []
loc = []
down = []
up = []
i = 0

Filt = filter()

F0 = 9
t2 = F0/10e6
Tout1 = t2

Yma = 0
buf = []
buf_sigm = []

bath_inf = bath()
sigm = 20
T1 = time.time()
while i <= 3000:
	t1 = 0
	# t1 = random.gauss(0,5e-8)
	# if i >= 1500 and i <= 2500:
	# 	t1 = random.gauss(0,100e-8)
	# if i >= 0 and i <= 2000:
	# 	F0 = F0 + 0.0001
	# if i >= 2000 and i <= 6000:
	# 	F0 = F0 - 0.0001
	# if i >= 6000 and i <= 8000:
	# 	F0 = F0 + 0.0001

	registers.pd_error =PHD(t1,t2)
	registers.gps_pll_lock  = lock_detector(bath_inf, registers.pd_error)

	registers.gps_dac_value = filterFanction(registers.pd_error, registers.gps_pll_lock , Filt)
	if (registers.pd_error >= 2*sigm or registers.pd_error <= -2*sigm) and registers.gps_pll_lock  == 1:
		registers.gps_dac_value = Yma
	
	if registers.gps_pll_lock:
		sigm = sigma(buf_sigm, registers.pd_error)
		Yma = moving_average(registers.gps_dac_value, buf)

	t2 = VCO(registers.gps_dac_value, F0, Tout1)
	Tout1 = t2

	out_err.append(t2)
	t.append(i)
	down.append(-1e-6)
	up.append(1e-6)
	loc.append(registers.gps_pll_lock *6*10e-7)
	input_err.append(t1)

	i = i + 1

Vvco = Yma
T2 = time.time()
T=T2-T1
print(T/3000)
ii = 0
hold_out = []
thold = []

t_last = out_err[len(out_err)-1]

# while ii < 50000:

# 	hold = VCO(Vvco, F0, t_last)
# 	t_last = hold
# 	hold_out.append(hold)
# 	thold.append(ii)
# 	ii = ii + 1

plt.plot(t,input_err)
plt.plot(t,out_err)
plt.plot(t,loc)
plt.plot(t,down)
plt.plot(t,up)
# plt.plot(tg,Gpfd)
# plt.plot(thold, hold_out)
plt.grid()
plt.show()
