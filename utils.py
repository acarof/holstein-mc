import numpy as np
import random


def write_file(list_, title):
    with open(title, 'w') as file_:
        file_.write('\n'.join(map(str, list_)) + '\n')


def write_traj(list_, title):
    with open(title, 'w') as file_:
        file_.write('\n'.join(['  '.join(map(str, _)) for _ in list_]) + '\n')


def create_hamiltonian(M, couplings, frequency, mass, xi, massy=0, freqy=0, beta=0, yi=[]):
    size = len(xi)
    if len(yi) == 0:
        yi = np.zeros(size)
    hamiltonian = couplings
    sum_ = 0.0
    sum_ = 0.5 * (mass * frequency ** 2) * np.sum(np.power(xi, 2)) + \
           0.5 * (massy * freqy ** 2) * np.sum(np.power(yi, 2))
    for index, x in enumerate(xi):
        hamiltonian[index, index] = M * x + sum_
        if (index + 1) < size:
            hamiltonian[index, index + 1] = beta * (yi[index + 1] - yi[index])
            hamiltonian[index + 1, index] = beta * (yi[index + 1] - yi[index])
    # print hamiltonian
    return hamiltonian


def give_eigen(hamiltonian, state):
    eigenValues, eigenVectors = np.linalg.eig(hamiltonian)
    idx = eigenValues.argsort()
    eigenValues = eigenValues[idx]
    eigenVectors = eigenVectors[:, idx]
    return eigenValues[state], eigenVectors[:, state]


def couplings_hamilt(size, dim, dim_coupling):
    network = np.zeros((size,) * dim)
    couplings = np.zeros((size ** dim, size ** dim))

    it = np.nditer(network, flags=['multi_index'])
    index = 0

    while not it.finished:
        index2 = 0
        it2 = np.nditer(network, flags=['multi_index'])
        while not it2.finished:
            vect = np.array(it2.multi_index) - np.array(it.multi_index)
            if (np.linalg.norm(vect)) == 1:
                pos = np.argmax(np.abs(vect))
                couplings[index, index2] = dim_coupling[pos]
            it2.iternext()
            index2 += 1
        it.iternext()
        index += 1
    return couplings


def run_monte_carlo(mass, dim_couplings, frequency, reorga,
                    massy, freqy, beta,
                    kbT, length, nsteps, grid, dim,
                    sampling=['local'],
                    toprint=['energies']):
    M = np.sqrt(reorga * mass * frequency ** 2)
    xunit = M / (mass * frequency ** 2)

    delta = xunit * grid
    properties = {}
    for prop in toprint:
        properties[prop] = []

    size = length ** dim
    couplings = couplings_hamilt(length, dim, dim_couplings)

    # print couplings

    def add_properties(properties, now_energy, now_vector, now_state, now_x):
        if 'energies' in properties:
            properties['energies'].append(now_energy)
        if 'IPR' in properties:
            properties['IPR'].append(1 / np.sum(np.power(now_vector, 4)))
        if 'state' in properties:
            properties['state'].append(now_state)
        if 'traj' in properties:
            properties['traj'].append(now_x)
        if 'trajy' in properties:
            properties['trajy'].append(now_y)

    now_x = np.random.normal(loc=0, scale=((mass * frequency ** 2) / kbT), size=size)
    now_y = np.random.normal(loc=0, scale=((mass * frequency ** 2) / kbT), size=size)
    now_x = np.zeros(size)
    now_state = 0
    hamiltonian = create_hamiltonian(M, couplings, frequency, mass, now_x, massy, freqy, beta, now_y)
    now_energy, now_vector = give_eigen(hamiltonian, now_state)
    add_properties(properties, now_energy, now_vector, now_state, now_x)

    for step in range(nsteps):
        do_test = True
        tosample = np.random.choice(sampling)
        try_state = now_state
        try_x = now_x
        try_y = now_y
        if tosample == 'local':
            try_x = now_x + np.random.normal(loc=0, scale=(delta), size=size)
        elif tosample == 'nonlocal':
            try_y = now_y + np.random.normal(loc=0, scale=(delta), size=size)
        elif tosample == 'state':
            try_state = now_state + np.random.choice([-1, 1])
            if try_state < 0 or try_state >= length:
                do_test = False
        if do_test:
            hamiltonian = create_hamiltonian(M, couplings, frequency, mass, try_x, massy, freqy, beta, try_y)
            try_energy, try_vector = give_eigen(hamiltonian, try_state)
            energy_diff = now_energy - try_energy
            boltmzan_weight = np.exp(energy_diff / kbT)
            if np.random.uniform() < boltmzan_weight:
                now_energy = try_energy
                now_vector = try_vector
                if tosample == 'local':
                    now_x = try_x
                elif tosample == 'state':
                    now_state = try_state
                elif tosample == 'nonlocal':
                    now_y = try_y

        add_properties(properties, now_energy, now_vector, now_state, now_x)
    return properties

