from flask import Flask
from flask import render_template
from flask import request
from flask import redirect, url_for
from numerize import numerize
import numpy as np

H = 69600/(3.086e22)

z = 0
omega_m = 0
omega_r = 0
omega_l = 0
omega_k = 0
c = 3e8

def hubble_func(k, o_m, o_r, o_l, o_k) -> float:
    first = o_m * (1+k)**3
    second = o_r * (1+k)**4
    third = o_l
    fourth = o_k * (1+k)**2

    overall = (1+k) * (first + second + third + fourth)**0.5

    return 1/overall

def co_rad(k, o_m, o_r, o_l, o_k) -> float:
    first = o_m * (1+k)**3
    second = o_r * (1+k)**4
    third = o_l
    fourth = o_k * (1+k)**2

    overall = (first + second + third + fourth)**0.5

    return 1/overall

def simpson1(z, o_m, o_r, o_l, o_k, h0, h) -> float:
    values = np.arange(0,z,h)
    sum = 0
    for i in range(len(values)):
        if i==0 or i==len(values)-1:
            sum += hubble_func(values[i], o_m, o_r, o_l, o_k)
        elif (i%3)==0:
            sum += 2 * hubble_func(values[i], o_m, o_r, o_l, o_k)
        else:
            sum += 3 * hubble_func(values[i], o_m, o_r, o_l, o_k)

    return (0.375 * h * sum)/(365 * 86400 * h0)

def simpson2(z, o_m, o_r, o_l, o_k, h0, h) -> float:
    values = np.arange(0,z,h)
    sum = 0
    for i in range(len(values)):
        if i==0 or i==len(values)-1:
            sum += co_rad(values[i], o_m, o_r, o_l, o_k)
        elif (i%3)==0:
            sum += 2 * co_rad(values[i], o_m, o_r, o_l, o_k)
        else:
            sum += 3 * co_rad(values[i], o_m, o_r, o_l, o_k)

    return (c * 0.375 * h * sum)/(3.086e16 * h0)

def ang_dist(z, o_m, o_r, o_l, o_k, h0, h) -> float:
    values = np.arange(0,z,h)
    sum = 0
    for i in range(len(values)):
        if i==0 or i==len(values)-1:
            sum += co_rad(values[i], o_m, o_r, o_l, o_k)
        elif (i%3)==0:
            sum += 2 * co_rad(values[i], o_m, o_r, o_l, o_k)
        else:
            sum += 3 * co_rad(values[i], o_m, o_r, o_l, o_k)

    return (c * 0.375 * h * sum)/( (1+z) * 3.086e16 * h0)

def lum_dist(z, o_m, o_r, o_l, o_k, h0, h) -> float:
    values = np.arange(0,z,h)
    sum = 0
    for i in range(len(values)):
        if i==0 or i==len(values)-1:
            sum += co_rad(values[i], o_m, o_r, o_l, o_k)
        elif (i%3)==0:
            sum += 2 * co_rad(values[i], o_m, o_r, o_l, o_k)
        else:
            sum += 3 * co_rad(values[i], o_m, o_r, o_l, o_k)

    return (c * (1+z) * 0.375 * h * sum)/(3.086e16 * h0)

def current_age(o_m, o_r, o_l, o_k, h0) -> float:

    value = np.arange(0,5000,0.1)
    sum = 0
    for i in range(len(value)):
        if i==0 or i==len(value)-1:
            sum += hubble_func(value[i], o_m, o_r, o_l, o_k)
        elif (i%3)==0:
            sum += 2 * hubble_func(value[i], o_m, o_r, o_l, o_k)
        else:
            sum += 3 * hubble_func(value[i], o_m, o_r, o_l, o_k)

    return (0.375 * 0.1 * sum)/(h0 * 365 * 86400)

def temp(z) -> float:
    return 2.73 * (1+z)

def ang_scale(z, o_m, o_r, o_l, o_k, h0, h) -> float:
    get = ang_dist(z, o_m, o_r, o_l, o_k, h0, h)
    circum = 2 * np.pi * get
    scale = circum/(360 * 3600)
    return scale

app = Flask(__name__, template_folder='template')

@app.route('/')
def hello():
    return render_template('second.html')

@app.route('/calc')
def submit():
    return render_template('first.html')

@app.route('/calc', methods=['POST'])
def calc():
    global z
    global omega_m
    global omega_r
    global omega_l
    global omega_k
    global H

    if request.method=='POST':
        z = request.form.get('Redshift', type=float)
        omega_m = request.form.get('Matter Density', type=float)
        omega_r = request.form.get('Radiation Density', type=float)
        omega_l = request.form.get('Dark Energy Density', type=float)
        omega_k = request.form.get('Curvature Density', type=float)
        H = request.form.get('Hubble Constant', type=float)

    H0 = (H * 1e3)/(3.086e22)

    h = 0
    if z<=50:
        h = 0.0001
    else:
        h = 0.01


    return render_template('first.html', ca=numerize.numerize(current_age(omega_m, omega_r, omega_l, omega_k, H0)), lt=numerize.numerize(simpson1(z, omega_m, omega_r, omega_l, omega_k, H0, h)) ,red=z, ra = numerize.numerize(current_age(omega_m, omega_r, omega_l, omega_k, H0) - simpson1(z, omega_m, omega_r, omega_l, omega_k, H0, h)), cd=numerize.numerize(simpson2(z, omega_m, omega_r, omega_l, omega_k, H0, h)), ad=numerize.numerize(ang_dist(z, omega_m, omega_r, omega_l, omega_k, H0, h)), ld=numerize.numerize(lum_dist(z, omega_m, omega_r, omega_l, omega_k, H0, h)), tz=temp(z), ags=numerize.numerize(ang_scale(z, omega_m, omega_r, omega_l, omega_k, H0, h)))

if __name__=='__main__':
    app.run(debug=True)

