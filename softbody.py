import numpy as np
import matplotlib.pyplot as plt
import os

path = 'frames'
os.system(f'mkdir {path}')

# parameters
N = 8                                       # grid size (N x N)
mass = 0.05                         
spring_k = 100.0                            # spring constant
spring_k_diagonal = spring_k / 2            # diagonal spring constant (smaller)
damping = 0.5                               # low damping = more jello
rest_length = 1.0                           
rest_length_diag = np.sqrt(2) * rest_length
time_step = 0.01                            # if unstable, try smaller dt 
num_steps = 1000                            # number of time steps
gravity = np.array([0, -9.81])
floor_y = 0                                 # floor y coordinate  
restitution = 0.75 
friction_coefficient = 100
yshift = 6                                  # sets height where body is dropped

# initialize positions
positions = np.array([(i, j) for i in range(N) for j in range(N)], dtype=np.float64)
positions[:, 1] += yshift
velocities = np.zeros_like(positions)

def find_neighbors(index):
    neighbors = {'ortho': [], 'diag': []}
    x, y = index % N, index // N
    for dx in [-1, 0, 1]:
        for dy in [-1, 0, 1]:
            if dx == 0 and dy == 0:
                continue  
            nx, ny = x + dx, y + dy
            if 0 <= nx < N and 0 <= ny < N:
                if dx == 0 or dy == 0:
                    neighbors['ortho'].append(ny * N + nx)
                else:
                    neighbors['diag'].append(ny * N + nx)
    return neighbors

def compute_spring_forces(positions):
    forces = np.zeros_like(positions)
    for i, pos in enumerate(positions):
        neighbors = find_neighbors(i)
        for j in neighbors['ortho']:
            direction = positions[j] - pos
            distance = np.linalg.norm(direction)
            if distance > 0:
                force_magnitude = spring_k * (distance - rest_length)
                force = (direction / distance) * force_magnitude
                forces[i] += force
                forces[j] -= force  # Newton's third law
        for j in neighbors['diag']:
            direction = positions[j] - pos
            distance = np.linalg.norm(direction)
            if distance > 0:
                force_magnitude = spring_k_diagonal * (distance - rest_length_diag)
                force = (direction / distance) * force_magnitude
                forces[i] += force
                forces[j] -= force
    return forces

def compute_damping_forces(positions, velocities):
    forces = np.zeros_like(positions)
    for i, pos in enumerate(positions):
        neighbors = find_neighbors(i)
        for j in neighbors['ortho'] + neighbors['diag']:
            relative_velocity = velocities[i] - velocities[j]
            damping_force = damping * relative_velocity
            forces[i] -= damping_force / len(neighbors['ortho'] + neighbors['diag'])
    return forces

def floor_collision(positions, velocities):
    for i, (x, y) in enumerate(positions):
        if y < floor_y:                                             # collision detection
            positions[i, 1] = floor_y                               # reset y position to floor level
            if velocities[i, 1] < 0:                                # only reverse velocity if it's moving towards the floor
                velocities[i, 1] = -velocities[i, 1] * restitution  # bounce effect
            
            normal_force = mass * abs(gravity[1])                   # normal force due to gravity
            friction_force_magnitude = friction_coefficient * normal_force
            velocity_magnitude = np.linalg.norm(velocities[i])
            
            if velocity_magnitude > 1e-5:
                friction_force = friction_force_magnitude * (velocities[i] / velocity_magnitude)
                velocities[i] -= friction_force * time_step

for step in range(num_steps):
    print(f'{step}/{num_steps} done')
    spring_forces = compute_spring_forces(positions)
    damping_forces = compute_damping_forces(positions, velocities)
    total_forces = spring_forces + damping_forces + mass * gravity
    velocities += (total_forces / mass) * time_step
    positions += velocities * time_step
    floor_collision(positions, velocities)
    #total_energy = calculate_energy(positions, velocities)
    #angular_momentum = calculate_angular_momentum(positions, velocities)
    
    fig, ax = plt.subplots()
    plt.scatter(positions[:, 0], positions[:, 1], color='blue', s=10)   # draw balls
    plt.plot([-16, N+15], [floor_y, floor_y], 'black')                  # draw floor
    plt.xlim(-16, N+15)
    plt.ylim(-16, N+15)
    plt.savefig(f'{path}/frame_{step:05}.png', bbox_inches='tight', pad_inches=0)
    plt.close("fig")
