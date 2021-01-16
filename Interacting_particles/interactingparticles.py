import math
import random 
import numpy as np
import matplotlib.pyplot as plt


# Function returns the distance of k particles to particle i 
# If the plains are mirrord round the middel plain 
def Smallest_distance(L_pos,k,i, box):
    distance = L_pos[k] - L_pos[i]
    if abs(distance) > box / 2 and L_pos[i] > 0:
        distance = box + distance
    if abs(distance) > box / 2 and L_pos[i] < 0:
        distance = distance - box
    return distance

## Special periodic conditions, making sure particle stays in the box
# Particle reapears on onther side box
def Pacman_boundery_conditions(nieuw, min, max, box):
        # Position is greater than max bos
        while nieuw < min:
            nieuw = nieuw + box
        # Position is smaller than min box
        while nieuw > max: 
            nieuw = nieuw - box 
        return nieuw

# Function returns the avarige energy over velocitys coordinates in L_v_vector_x and L_v_vector_y  
def E_kin_avr(L_v_vector_x, L_v_vector_y):
    E_kin_tot = 0
    for i in range(0, len(L_v_vector_x)):
        E_kin = 0.5 * ((L_v_vector_x[i])**2 + (L_v_vector_y[i])**2)
        E_kin_tot += E_kin
    
    E_kin_avr = E_kin_tot / len(L_v_vector_x)
    return E_kin_avr

# Calculating potential energy of particle i due to particle k
def Potential(r_k, V_0, r_0):
    V_k = 4 * V_0 * ((r_0 / r_k)**12 - (r_0 / r_k)**6)
    return V_k

# Calculating avaerige potential energy of i particles 
def E_pot_avr(L_V_i):
    Pot_tot = 0
    for i in range(0,len(L_V_i)):
        Pot_tot += L_V_i[i]
    Pot_avr = Pot_tot / len(L_V_i)
    return Pot_avr

# Function returns the x,y components of the acceleration vector for particle i due to particle k 
# Imput is: r and angle components of the distance between i and k particle. T
# the V_o and r_0 of the Lennard_Jones potential
# And the mass of particle i
def Acceleration_k(r_k, angle_k, V_0, r_0, mass_i):
    a_x = (1/mass_i) * (V_0 * (24 / r_k) * (2 * (r_0 / r_k)**12 - (r_0 / r_k)**6)) * math.cos(angle_k)
    a_y = (1/mass_i) * (V_0 * (24 / r_k) * (2 * (r_0 / r_k)**12 - (r_0 / r_k)**6)) * math.sin(angle_k)
    return [a_x, a_y]

# Function wil return the x,y coordinates of the net acceleration on particle i and the potential of particle i
# due to the Lennard_Jones potential force of al k particles 
# Function uses prviously defined functions: Acceleration_k, Smallest_distance and Potential
def Net_acceleration_i(L_pos_x, L_pos_y, x_box, y_box, i, V_0, r_0, mass_i):
    net_a_x = 0
    net_a_y = 0

    # Defining parameter for potential of particle i
    V_i = 0

    for k in range(0, len(L_pos_x)):
        # Checking force due to al particles k (k is not = to i)
        if k != i:
            # First checking distances on the plain
            d_x_k = Smallest_distance(L_pos_x,k,i, x_box)
            d_y_k = Smallest_distance(L_pos_y,k,i, y_box)

            # Determining radial distance 
            r_k = math.sqrt(d_x_k**2 + d_y_k**2)
            
            # Calculating potential due to particle k
            V_k = Potential(r_k, V_0, r_0)
            V_i += V_k

            # Determining angle to x ax
            angle_k = np.arctan2(d_y_k, d_x_k)
            if angle_k < 0:
                angle_k = 2 * math.pi + angle_k

            net_a_x += Acceleration_k(r_k, angle_k, V_0, r_0, mass_i)[0]
            net_a_y += Acceleration_k(r_k, angle_k, V_0, r_0, mass_i)[1]

    return[net_a_x, net_a_y, V_i]
         
## Function returns list of starting positions, velocitys 
# and calculating the starting accelerations and potentials
# Function uses prviously defined functions: Net_acceleration_i, 
def Start_values(lowbound, upbound, step, v_max, x_box, y_box,  V_0, r_0, mass_i):

    L_pos_x = []
    L_pos_y = []
    L_v_vector_x = []
    L_v_vector_y = []
    L_a_vector_x = []
    L_a_vector_y = []
    
    L_V_i = []

    # Position and velocity
    for pos_x in np.arange(lowbound,upbound + step, step):
        for pos_y in np.arange(lowbound,upbound + step, step):
            L_pos_x.append(pos_x)
            L_pos_y.append(pos_y)

    # Adding random displacements
    for i in range(0, len(L_pos_x)):
        # Adding random displacements
        random_displacement = r_0 / 3
        L_pos_x[i] = L_pos_x[i] + 2 * (random.random() - 0.5) * random_displacement
        L_pos_y[i] = L_pos_y[i] + 2 * (random.random() - 0.5) * random_displacement

        ##  Making random staring velocitys 
        # Making random speed angle 
        angle = (2 * math.pi) * random.random()
        # Making random speed between 0 and v_max
        speed = v_max * random.random()
        # Making list of speed vector components
        v_vector_x = speed * math.cos(angle)
        v_vector_y = speed * math.sin(angle)

        L_v_vector_x.append(v_vector_x)
        L_v_vector_y.append(v_vector_y)

        ## Determining starting acceleration (# Zou 0 moeten zijn)
        net_a_x = Net_acceleration_i(L_pos_x, L_pos_y, x_box, y_box, i, V_0, r_0, mass_i)[0]
        net_a_y = Net_acceleration_i(L_pos_x, L_pos_y, x_box, y_box, i, V_0, r_0, mass_i)[1]

        L_a_vector_x.append(net_a_x)
        L_a_vector_y.append(net_a_y)

        ## Making list for potential energy of every particle 
        V_i = Net_acceleration_i(L_pos_x, L_pos_y, x_box, y_box, i, V_0, r_0, mass_i)[2]
        L_V_i.append(V_i)

    # Returning starting values
    return [L_pos_x, L_pos_y, L_v_vector_x, L_v_vector_y, L_a_vector_x, L_a_vector_y, L_V_i]
            
# Function returns lists with new values for: L_pos_x, L_pos_y, L_v_vector_x, L_v_vector_y, 
# L_a_vector_x, L_a_vector_y and L_V_i after t_step using Verlet algorithm
# Function usus pryviosly defind functions: Pecman_boundery_conditions and Net_acceleration_i
def New_xy_v_a(L_pos_x, L_pos_y,L_v_vector_x,L_v_vector_y,L_a_vector_x, L_a_vector_y,L_V_i,t_step,x_box, y_box, V_0, r_0, mass_i):

    # Coordinates for the box
    x_min = - x_box / 2
    x_max = x_box / 2
    y_min =  - y_box / 2
    y_max = y_box / 2

    # Calculating new position, velocity and acceleration for every particle i
    for i in range(0, len(L_pos_x)):
        x_nieuw = L_pos_x[i] + L_v_vector_x[i] * t_step + 0.5 * L_a_vector_x[i] * (t_step)**2
        y_nieuw = L_pos_y[i] + L_v_vector_y[i] * t_step + 0.5 * L_a_vector_y[i] * (t_step)**2

        x_nieuw = Pacman_boundery_conditions(x_nieuw, x_min, x_max, x_box)
        y_nieuw = Pacman_boundery_conditions(y_nieuw, y_min, y_max, y_box)

        L_pos_x[i] = x_nieuw
        L_pos_y[i] = y_nieuw

        a_net_x_old = L_a_vector_x[i]
        a_net_y_old = L_a_vector_y[i]

        a_net_x_nieuw = Net_acceleration_i(L_pos_x, L_pos_y, x_box, y_box, i, V_0, r_0, mass_i)[0]
        a_net_y_nieuw = Net_acceleration_i(L_pos_x, L_pos_y, x_box, y_box, i, V_0, r_0, mass_i)[1]

        L_a_vector_x[i] = a_net_x_nieuw
        L_a_vector_y[i] = a_net_y_nieuw

        v_x_nieuw = L_v_vector_x[i] + 0.5 * (a_net_x_nieuw + a_net_x_old)
        v_y_nieuw = L_v_vector_y[i] + 0.5 * (a_net_y_nieuw + a_net_y_old)

        L_v_vector_x[i] = v_x_nieuw
        L_v_vector_y[i] = v_y_nieuw

        L_V_i[i] = Net_acceleration_i(L_pos_x, L_pos_y, x_box, y_box, i, V_0, r_0, mass_i)[2]
    
    return[L_pos_x, L_pos_y,L_v_vector_x,L_v_vector_y,L_a_vector_x, L_a_vector_y]

def Simulation(t_max, t_step, v_max_start, x_box, y_box, V_0, r_0, mass_i):
    
    # Starting by making Lists of starting values where 16 particles are on a even grid
    lowbound = -1.5
    upbound = 1.5
    step = 1

    # Coordinates for the box
    x_min = - x_box / 2
    x_max = x_box / 2
    y_min =  - y_box / 2
    y_max = y_box / 2


    #L_pos_x, L_pos_y, L_v_vector_x, L_v_vector_y, L_a_vector_x, L_a_vector_y
    L_pos_x = Start_values(lowbound, upbound, step, v_max_start, x_box, y_box, V_0, r_0, mass_i)[0]
    L_pos_y = Start_values(lowbound, upbound, step, v_max_start, x_box, y_box, V_0, r_0, mass_i)[1]

    
    L_v_vector_x = Start_values(lowbound, upbound, step, v_max_start, x_box, y_box, V_0, r_0, mass_i)[2]
    L_v_vector_y = Start_values(lowbound, upbound, step, v_max_start, x_box, y_box, V_0, r_0, mass_i)[3]

    L_a_vector_x = Start_values(lowbound, upbound, step, v_max_start, x_box, y_box, V_0, r_0, mass_i)[4]
    L_a_vector_y = Start_values(lowbound, upbound, step, v_max_start, x_box, y_box, V_0, r_0, mass_i)[5]

    L_V_i = Start_values(lowbound, upbound, step, v_max_start, x_box, y_box, V_0, r_0, mass_i)[6]

    # Defining listst
    L_time = []
    L_E_kin_avarige = []
    L_E_pot_avarige = []

    for t in np.arange(0, t_max, t_step):

        # Determining avaeruge kinetic energy
        E_kin_avarige = E_kin_avr(L_v_vector_x, L_v_vector_y)
        E_pot_avarige = E_pot_avr(L_V_i)
        L_E_kin_avarige.append(E_kin_avarige)
        L_E_pot_avarige.append(E_pot_avarige)
        L_time.append(t)

        # Determinung new position, velocity and acceleration for t+t_step
        # Retturn: L_pos_x, L_pos_y,L_v_vector_x,L_v_vector_y,L_a_vector_x, L_a_vector_y
        L_pos_x =  New_xy_v_a(L_pos_x, L_pos_y,L_v_vector_x,L_v_vector_y,L_a_vector_x, L_a_vector_y,L_V_i,t_step,x_box, y_box, V_0, r_0, mass_i)[0]
        L_pos_y =  New_xy_v_a(L_pos_x, L_pos_y,L_v_vector_x,L_v_vector_y,L_a_vector_x, L_a_vector_y,L_V_i,t_step,x_box, y_box, V_0, r_0, mass_i)[1]

        L_v_vector_x = New_xy_v_a(L_pos_x, L_pos_y,L_v_vector_x,L_v_vector_y,L_a_vector_x, L_a_vector_y,L_V_i,t_step,x_box, y_box, V_0, r_0, mass_i)[2]
        L_v_vector_y = New_xy_v_a(L_pos_x, L_pos_y,L_v_vector_x,L_v_vector_y,L_a_vector_x, L_a_vector_y,L_V_i,t_step,x_box, y_box, V_0, r_0, mass_i)[3]

        L_a_vector_x = New_xy_v_a(L_pos_x, L_pos_y,L_v_vector_x,L_v_vector_y,L_a_vector_x, L_a_vector_y,L_V_i,t_step,x_box, y_box, V_0, r_0, mass_i)[4]
        L_a_vector_y = New_xy_v_a(L_pos_x, L_pos_y,L_v_vector_x,L_v_vector_y,L_a_vector_x, L_a_vector_y,L_V_i,t_step,x_box, y_box, V_0, r_0, mass_i)[5]


        # Plotting position of i particles at time t
        plt.plot(L_pos_x, L_pos_y, 'ro')
        # Plotting acceleration vector of i particles at time t
        plt.quiver(L_pos_x,L_pos_y,L_a_vector_x,L_a_vector_y,angles='xy', scale_units='xy', scale=1)
        
        # Plotting box
        plt.xlim(x_min,x_max)
        plt.ylim(y_min,y_max)

        # Plotttin time in title 
        plt.text(x_min + 0.5, y_min + 0.5, "({})".format(t))

        # Drawing plot 
        plt.draw()
        plt.pause(0.01)
        plt.clf()
    
    plt.plot(L_time,L_E_kin_avarige)
    plt.plot(L_time,L_E_pot_avarige)
    plt.ylabel("Kinetic energy [J]")
    plt.xlabel("time [s]")
    plt.savefig("Avarige kinetic energy")

Simulation(0.004,0.00001,0.1,4,4,1,1,1)

### Things that go wrong: 
# 1. E_kin = 0 at t = 0, this is not correct becaus v => 0 at t = 0
# 2. V_pot = 0 for al t
# 3. If step between particle starting positions = r_0 (Question 2) the system dous not work


## Queestions:
# What is a good timescale 
# Are the boundery conditions applyd correctly in Pecman_boundery_conditions? 
# Is the verlet algorithm applyd correctly in New_xy_v_a?

## To Do: 
# Speed and velocity distributions