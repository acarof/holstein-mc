import numpy as np


def write_file(list_, title):
    with open(title, 'w') as file_:
        file_.write('\n'.join(map(str, list_)))


def create_hamiltonian(M, coupling, frequency, mass, *xi):
    size = len(xi)
    hamiltonian = np.zeros([size, size])
    sum_ = 0.0
    for x in xi:
        sum_ += 0.5 * (mass * frequency ** 2) * x ** 2
    for index, x in enumerate(xi):
        hamiltonian[index, index] = M * x + sum_
        if (index + 1) < size:
            hamiltonian[index, index + 1] = coupling
            hamiltonian[index + 1, index] = coupling
    # print hamiltonian
    return hamiltonian




def give_lowest_eigen(hamiltonian):
    eigenValues, eigenVectors = np.linalg.eig(hamiltonian)
    idx = eigenValues.argsort()[::-1]
    eigenValues = eigenValues[idx]
    eigenVectors = eigenVectors[:, idx]
    return eigenValues[0], eigenVectors[0]


def run_monte_carlo(mass, coupling, frequency, reorga, kbT, size, nsteps, grid):
    M = np.sqrt(reorga * mass * frequency ** 2 / 2)
    xunit = M / (mass * frequency ** 2)

    delta = xunit / grid

    now_x = np.random.normal(loc=0, scale=((mass * frequency ** 2) / kbT), size=size)
    hamiltonian = create_hamiltonian(M, coupling, frequency, mass, *now_x)
    now_energy, now_vector = give_lowest_eigen(hamiltonian)
    now_ipr = 1 / np.sum(np.power(now_vector, 4))
    energies = [now_energy]
    ipr = [now_ipr]

    for step in range(nsteps):
        try_x = now_x + np.random.normal(loc=0, scale=(delta), size=size)
        hamiltonian = create_hamiltonian(M, coupling, frequency, mass, *try_x)
        try_energy, try_vector = give_lowest_eigen(hamiltonian)
        try_ipr = 1 / np.sum(np.power(try_vector, 4))

        energy_diff = energies[-1] - try_energy
        boltmzan_weight = np.exp(energy_diff / kbT)
        random = np.random.uniform()
        if random < boltmzan_weight:
            now_energy = try_energy
            now_x = try_x
            now_ipr = try_ipr
        else:
            now_energy = energies[-1]
            now_x = now_x
            now_ipr = now_ipr
        energies.append(now_energy)
        ipr.append(now_ipr)
    return energies, ipr

