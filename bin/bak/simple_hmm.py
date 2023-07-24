#! /usr/bin/env python
from __future__ import print_function
import numpy as np
from cStringIO import StringIO
import h5py
import os, sys

def prepare_output_file(path):
    if path.startswith('/dev'):
        return
    dirpath = os.path.dirname(path)
    if not os.path.exists(dirpath):
        os.makedirs(dirpath)

class GaussianHMM(object):
    @staticmethod
    def load_model(filename):
        model = GaussianHMM()
        model.load_weights(filename)
        return model

    def __init__(self, n_states=2):
        self.n_states = n_states
        self.p_state = np.zeros(self.n_states)
        self.p_transition = np.zeros((self.n_states, self.n_states))
        self.p_emission_mean = np.zeros(self.n_states)
        self.p_emission_std = np.zeros(self.n_states)
        self.init_weights()

    def init_weights(self):
        from scipy.stats import norm
        self.p_state = np.random.uniform(size=self.n_states)
        self.p_state = np.exp(self.p_state*10.0)
        self.p_state /= self.p_state.sum()
        self.p_transition = np.random.uniform(size=(self.n_states, self.n_states))
        self.p_transition = np.exp(self.p_transition*10.0)
        self.p_transition /= self.p_transition.sum(axis=1).reshape((-1, 1))
        self.p_emission_mean = np.random.uniform(-10, 10, size=self.n_states)
        self.p_emission_std = np.ones(self.n_states)
        self.emission = norm

    def init_weights_kmeans(self, x, pseudo_counts_state=1, pseudo_counts_transition=1):
        from sklearn.cluster import KMeans
        # remove missing values
        x = x[np.logical_not(np.isnan(x))]
        kmeans = KMeans(self.n_states)
        kmeans.fit(x.reshape((-1, 1)))
        self.p_emission_mean = np.zeros(self.n_states)
        self.p_emission_std = np.zeros(self.n_states)
        self.p_state = np.zeros(self.n_states)
        for k in range(self.n_states):
            x_cluster = x[kmeans.labels_ == k]
            self.p_state[k] = len(x_cluster)
            self.p_emission_mean[k] = np.mean(x_cluster)
            self.p_emission_std[k] = np.std(x_cluster)
        self.p_state += pseudo_counts_state
        self.p_state /= np.sum(self.p_state)

        transition_counts = np.full((self.n_states, self.n_states), pseudo_counts_transition, dtype='float64')
        for t in range(1, len(x)):
            transition_counts[kmeans.labels_[t - 1], kmeans.labels_[t]] += 1
        self.p_transition = transition_counts / transition_counts.sum(axis=1).reshape((-1, 1))

        print('Initial state probabilities = ', self.p_state)
        print('Initial transition matrix = ', self.p_transition)
        print('Initial emission mean = ', self.p_emission_mean)
        print('Initial emission std = ', self.p_emission_std)

    def save(self, filename):
        self.save_weights(filename)

    def save_weights(self, filename):
        f = h5py.File(filename, 'w')
        f.create_dataset('p_state', data=self.p_state)
        f.create_dataset('p_transition', data=self.p_transition)
        f.create_dataset('p_emission_mean', data=self.p_emission_mean)
        f.create_dataset('p_emission_std', data=self.p_emission_std)
        f.close()

    def load_weights(self, filename):
        f = h5py.File(filename, 'r')
        self.p_state = f['p_state'][:]
        self.n_states = self.p_state.shape[0]
        self.p_emission_mean = f['p_emission_mean'][:]
        self.p_emission_std = f['p_emission_std'][:]
        self.p_transition = f['p_transition'][:]
        f.close()

    def __str__(self):
        f = StringIO()
        print('Initial state probabilities = ', self.p_state, file=f)
        print('Transition probabilities = ', self.p_transition, file=f)
        print('Emission mean = ', self.p_emission_mean, file=f)
        print('Emission std = ', self.p_emission_std, file=f)
        return f.getvalue()

    def fit(self, x, max_iter=50, verbose=False):
        self.init_weights_kmeans(x)
        likelihood = 0
        for i in range(max_iter):
            gamma, xi, new_likelihood = self.estep(x)
            self.mstep(x, gamma, xi)
            if verbose:
                print('Iteration %d: log likelihood = %f'%(i + 1, likelihood))
            #if (new_likelihood > likelihood) or (likelihood == 0):
            if True:
                likelihood = new_likelihood
            else:
                break

    def likelihood(self, x):
        alpha_, beta_, c, e, p_x = self.forward_backward(x)
        return p_x

    def forward_backward(self, x):
        T = len(x)
        K = self.n_states
        A = self.p_transition
        alpha_ = np.zeros((T, K))
        beta_ = np.zeros((T, K))
        c = np.zeros(T)
        e = np.empty((T, K))
        pi = self.p_state

        # calculate emission probabilities
        nan_ind = np.nonzero(np.isnan(x))[0]
        x[nan_ind] = 0
        for k in range(K):
            e[:, k] = self.emission.pdf(x, self.p_emission_mean[k], self.p_emission_std[k])
        # ignore emission probabilities for missing values
        e[nan_ind, k] = 1.0

        # forward pass to calculate alpha
        for k in range(K):
            alpha_[0, k] = e[0, k]*pi[k]
        c[0] = np.sum(alpha_[0])
        alpha_[0] /= c[0]
        for t in range(1, T):
            for k in range(K):
                alpha_[t, k] = e[t, k]*np.sum(alpha_[t - 1]*A[:, k])
            c[t] = np.sum(alpha_[t])
            alpha_[t] /= c[t]

        # backward pass to calculate beta
        beta_[T - 1] = 1.0
        for t in range(T - 2, -1, -1):
            for k in range(K):
                beta_[t, k] = np.sum(beta_[t + 1]*e[t + 1]*A[k])
            beta_[t] /= c[t + 1]

        likelihood = np.sum(np.log(c))

        #print('Alpha = ', alpha_)
        #print('Beta = ', beta_)
        #print('Expected states = ', alpha_ * beta_)
        #print('Scale factor = ', c)
        #print('Likelihood = ', likelihood)

        return alpha_, beta_, c, e, likelihood

    def estep(self, x):
        T = len(x)
        K = self.n_states
        alpha_, beta_, c, e, likelihood = self.forward_backward(x)
        A = self.p_transition
        gamma = alpha_*beta_
        xi = np.empty((T - 1, K, K))
        for t in range(1, T):
            xi[t - 1] = 1.0/c[t]*(alpha_[t - 1].reshape((K, 1)))*(e[t].reshape((1, K)))*A*(beta_[t].reshape((1, K)))
        #print('Expected states: ', gamma)
        #print('Expected transitions: ', xi)
        return gamma, xi, likelihood

    def mstep(self, x, expected_states, expected_transitions):
        gamma = expected_states
        xi = expected_transitions
        K = self.n_states

        pi_new = np.zeros(K)
        A_new = np.zeros((K, K))
        mean_new = np.zeros(K)
        std_new = np.zeros(K)

        for k in range(K):
            A_new = np.sum(xi, axis=0)/np.sum(xi, axis=(0, 2)).reshape((K, 1))
            mean_new[k] = np.sum(gamma[:, k]*x)/np.sum(gamma[:, k])
            std_new[k] = np.sqrt(np.sum(gamma[:, k]*((x - mean_new[k])**2))/np.sum(gamma[:, k]))
        #print('New transition matrix = ', A_new)
        #print('New mean = ', mean_new)
        #print('New std = ', std_new)
        self.p_emission_mean = mean_new
        self.p_emission_std = std_new
        self.p_transition = A_new
        return mean_new, std_new

    def viterbi_decode(self, x):
        T = len(x)
        K = self.n_states
        pi = self.p_state
        A = self.p_transition
        omega = np.zeros((T, K))
        e = np.zeros((T, K))
        # calculate emission probabilities
        # ignore emission probabilities for missing values
        nan_ind = np.nonzero(np.isnan(x))[0]
        x[nan_ind] = 0
        for k in range(K):
            e[:, k] = self.emission.pdf(x, self.p_emission_mean[k], self.p_emission_std[k])
        e[nan_ind, k] = 1.0

        # recursion
        omega[0] = np.log(pi) + np.log(e[0])
        max_origin = np.zeros((T, K), dtype='int32')
        max_origin[0] = 0
        for t in range(1, T):
            val = np.log(A) + omega[t - 1].reshape((K, 1))
            max_origin[t] = np.argmax(val, axis=0)
            omega[t] = np.log(e[t]) + val[max_origin[t], np.arange(K)]
        # backtrack
        path = np.zeros(T, dtype='int32')
        path[T - 1] = np.argmax(omega[T - 1])
        for t in range(T - 2, -1, -1):
            path[t] = max_origin[t + 1, path[t + 1]]
        return path, omega[T - 1, path[T - 1]]

    def mea_decode(self, x):
        gamma, xi, likelihood = self.estep(x)
        return np.argmax(gamma, axis=1), np.max(gamma, axis=1)

    def sample(self, length):
        x = np.empty(length)
        z = np.empty(length, dtype='int32')
        z[0] = np.random.choice(self.n_states, p=self.p_state)
        for t in range(1, length):
            z[t] = np.random.choice(self.n_states, p=self.p_transition[z[t - 1]])
        for t in range(length):
            x[t] = self.emission.rvs(self.p_emission_mean[z[t]], self.p_emission_std[z[t]])
        return x, z

if __name__ == '__main__':
    np.random.seed(1235)
    generator = GaussianHMM(5)
    print('generator = ', generator)
    #print('x = [', ' '.join(['%.2f'%a for a in x]), ']')
    #print('z = ', z)
    model_file = 'trained_models/hmm/GaussianHMM'
    if os.path.exists(model_file):
        model = GaussianHMM.load_model(model_file)
    else:
        x, z = generator.sample(5000)
        x[np.random.choice(x.shape[0], size=int(x.shape[0]*0.1), replace=True)] = np.nan
        model = GaussianHMM(5)
        model.fit(x, max_iter=50, verbose=True)
        prepare_output_file(model_file)
        model.save(model_file)
    state_rank_generator = np.argsort(np.argsort(generator.p_emission_mean))
    state_rank_model = np.argsort(np.argsort(model.p_emission_mean))
    print('Generator state ranks: ', state_rank_generator)
    print('Model state ranks: ', state_rank_model)
    print('generator = ', generator)
    x, z = generator.sample(50)
    #x[np.random.choice(x.shape[0], size=5, replace=True)] = np.nan
    z = np.take(state_rank_generator, z)
    print('Generate new x = ', x)
    print('z = ', z)
    z_pred, likelihood = model.viterbi_decode(x)
    z_pred = np.take(state_rank_model, z_pred)
    print('Viterbi decoded states = ', z_pred)
    print('Viterbi likelihood = ', likelihood)
    print('Categorical accuracy = ', float(np.count_nonzero(z_pred == z))/len(x))
    z_pred, likelihood = model.mea_decode(x)
    z_pred = np.take(state_rank_model, z_pred)
    print('MEA decoded states = ', z_pred)
    print('MEA likelihood = [', ' '.join(['%.2f'%p for p in likelihood]), ']')