#!/usr/bin/env python3
from math import floor
import random
import pygame as pg

print("this is fluidoh :)")

N = 64  #size
iteration = 1
SCALE = 5
t = 0

def IX(x, y):   #return the index in the array
    x = constrain(x, 0, N-1)
    y = constrain(y, 0, N-1)
    return x + (y*N)

def set_bnd(b, x):
    for i in range(1, N-1):
        x[IX(i, 0  )] = -x[IX(i, 1  )] if b == 2 else x[IX(i, 1 )]
        x[IX(i, N-1)] = -x[IX(i, N-2)] if b == 2 else x[IX(i, N-2)]

    for j in range(1, N-1):
        x[IX(0, j)] = -x[IX(1, j)] if b == 1 else x[IX(1, j)]
        x[IX(N-1, j)] = -x[IX(N-2, j)] if b == 1 else x[IX(N-2, j)]

    x[IX(0, 0)] = 0.5 * (x[IX(1, 0)] + x[IX(0, 1)])
    x[IX(0, N-1)] = 0.5 * (x[IX(1, N-1)] + x[IX(0, N-2)])
    x[IX(N-1, 0)] = 0.5 * (x[IX(N-2, 0)] + x[IX(N-1, 1)])
    x[IX(N-1, N-1)] = 0.5 * (x[IX(N-2, N-1)] + x[IX(N-1, N-2)])

    return x

def lin_solve(b, x, x0, a, c):
    cRecip = 1.0 / c
    for k in range(iteration):
        for j in range(1, N-1):
            for i in range(1, N-1):
                x[IX(i, j)] = (  x0[IX(i, j)] + a*( x[IX(i+1, j)] + x[IX(i-1, j)] + x[IX(i, j+1)] + x[IX(i, j-1)]  ) ) * cRecip

        x = set_bnd(b, x)
    return x, x0

def diffuse(b, x, x0, diff, dt):
    a = dt * diff * (N - 2) * (N - 2)
    x, x0 = lin_solve(b, x, x0, a, 1 + 4 * a)
    return x, x0


def project(velocX, velocY, p, div):
    for j in range(1, N-1):
        for i in range(1, N-1):
            div[IX(i, j)] = -0.5*(velocX[IX(i+1, j)]-velocX[IX(i-1, j)]+velocY[IX(i, j+1)]-velocY[IX(i, j-1)])/N
            p[IX(i, j)] = 0   

    div = set_bnd(0, div) 
    p = set_bnd(0, p)
    p, div = lin_solve(0, p, div, 1, 4)

    for j in range(1, N-1):
        for i in range(1, N-1):
            velocX[IX(i, j)] -= 0.5 * (  p[IX(i+1, j)] - p[IX(i-1, j)]) * N
            velocY[IX(i, j)] -= 0.5 * (  p[IX(i, j+1)] - p[IX(i, j-1)]) * N

    set_bnd(1, velocX)
    set_bnd(2, velocY)
    return velocX, velocY, p, div

def advect(b, d, d0, velocX, velocY, dt):
  #float i0, i1, j0, j1;

    dtx = dt * (N - 2)
    dty = dt * (N - 2)

  #float s0, s1, t0, t1;
  #float tmp1, tmp2, x, y;

    Nfloat = N
  #float ifloat, jfloat
  #int i, j
    jfloat = 1.0
    for j in range(1, N-1):
        ifloat = 1.0
        for i in range(1, N-1):
            tmp1 = dtx * velocX[IX(i, j)]
            tmp2 = dty * velocY[IX(i, j)]
            x    = ifloat - tmp1 
            y    = jfloat - tmp2

            if (x < 0.5):
                x = 0.5 
            if (x > Nfloat + 0.5):
                x = Nfloat + 0.5 
            i0 = floor(x) 
            i1 = i0 + 1.0
            if (y < 0.5):
                y = 0.5 
            if (y > Nfloat + 0.5):
                y = Nfloat + 0.5 
            j0 = floor(y)
            j1 = j0 + 1.0 

            s1 = x - i0 
            s0 = 1.0 - s1 
            t1 = y - j0 
            t0 = 1.0 - t1

            i0i = int(i0)
            i1i = int(i1)
            j0i = int(j0)
            j1i = int(j1)

            #DOUBLE CHECK THIS!!!
            d[IX(i, j)] = s0 * (t0 * d0[IX(i0i, j0i)] + t1 * d0[IX(i0i, j1i)]) + s1 * (t0 * d0[IX(i1i, j0i)] + t1 * d0[IX(i1i, j1i)])

            ifloat += 1
        jfloat += 1

    d = set_bnd(b, d)

    return d, d0, velocX, velocY

def constrain(val, min_val, max_val):
    return min(max_val, max(min_val, val))


class Fluid:
    """
    int size;
    float dt;
    float diff;
    float visc;
    
    float[] s;
    float[] density;
    
    float[] Vx;
    float[] Vy;
    float[] Vz;

    float[] Vx0;
    float[] Vy0;
    float[] Vz0;
    """

    s = []
    density = []
    
    Vx = []
    Vy = []
    Vz = []

    Vx0 = []
    Vy0 = []
    Vz0 = []

    def __init__(self, dt, diffusion, viscosity):
        self.size = N
        self.dt = dt
        self.diff = diffusion
        self.visc = viscosity

        self.s = [0] * N * N
        self.density = [0] * N * N
        
        self.Vx = [0] * N * N
        self.Vy = [0] * N * N

        self.Vx0 = [0] * N * N
        self.Vy0 = [0] * N * N

    def addDensity(self, x, y, amount):
        index = IX(x, y)
        self.density[index] += amount
    
    def addVelocity(self, x, y, amountX, amountY):
        index = IX(x, y)
        self.Vx[index] += amountX
        self.Vy[index] += amountY

    def step(self):
        """
        N = self.size
        visc    = self.visc
        diff    = self.diff
        dt      = self.dt
        Vx      = self.Vx
        Vy      = self.Vy
        Vx0     = self.Vx0
        Vy0     = self.Vy0
        s       = self.s
        density = self.density
        """

        self.Vx0, self.Vx = diffuse(1, self.Vx0, self.Vx, self.visc, self.dt)
        self.Vy0, self.Vy = diffuse(2, self.Vy0, self.Vy, self.visc, self.dt)

        self.Vx0, self.Vy0, self.Vx, self.Vy = project(self.Vx0, self.Vy0, self.Vx, self.Vy)

        self.Vx, self.Vx0, self.Vx0, self.Vy0 = advect(1, self.Vx, self.Vx0, self.Vx0, self.Vy0, self.dt)
        self.Vy, self.Vy0, self.Vx0, self.Vy0 = advect(2, self.Vy, self.Vy0, self.Vx0, self.Vy0, self.dt)

        self.Vx, self.Vy, self.Vx0, self.Vy0 = project(self.Vx, self.Vy, self.Vx0, self.Vy0)

        self.s, self.density = diffuse(0, self.s, self.density, self.diff, self.dt)
        self.density, self.s, self.Vx, self.Vy = advect(0, self.density, self.s, self.Vx, self.Vy, self.dt)

    def renderD(self):
        #colorMode(HSB, 255);

        for i in range(N):
            for j in range(N):
                x = i * SCALE
                y = j * SCALE
                d = self.density[IX(i, j)]
                #fill((d + 50) % 255,200,d)
                #noStroke()
                #square(x, y, SCALE)
                pg.draw.rect(scr, (constrain(d+50, 0, 255), constrain(d+50, 0, 255), constrain(d+50, 0, 255)), pg.Rect(x, y, SCALE, SCALE))
    
    #add renderV and fadeD
    def fadeD(self):
        for i in range(0, len(self.density)):
            d = self.density[i]
            self.density[i] = constrain(d-0.02, 0, 255)


fluid = Fluid(1, 0, 0.0000001)
pg.init()
scr = pg.display.set_mode((N*SCALE, N*SCALE))
scr.fill((0,0,0))
running = True
while running: 
    for event in pg.event.get(): 
        if event.type == pg.QUIT: 
            running = False
    #width, height = pg.display.get_surface().get_size()
    cx = int(0.5*N)
    cy = int(0.5*N)

    for i in range(-1, 2):
        for j in range(-1, 2):
            fluid.addDensity(cx+i, cy+j, random.randrange(50, 150))

    #for i in range(0, 2):
        #angle = noise
    vx = 3 #random.randrange(-1, 1)
    vy = 0 #random.randrange(-1, 1)
    t += 0.01
    fluid.addVelocity(cx, cy, vx, vy)
    
    fluid.step()
    fluid.renderD()
    fluid.fadeD()


    pg.display.flip()
pg.quit()