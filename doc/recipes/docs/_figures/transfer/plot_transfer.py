#!/usr/bin/env python3


'''
Generation of activation function plots.

The following functions are visualized:
tanh, sigmoid, heaviside, relu, linear, gaussian
'''


import numpy as np
from matplotlib import pyplot as plt


def main():
    '''Main driver routine.'''

    plot_tanh()
    plot_sigmoid()
    plot_gaussian()
    plot_linear()
    plot_relu()
    plot_heaviside()


def plot_tanh():
    '''Plots the hyperbolic tangent in given interval.'''

    xx = np.linspace(-np.pi, np.pi, 1000)
    tanh = np.tanh(xx)

    plt.figure(1, figsize=[7, 5])
    plt.title('Hyperbolic Tangent')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$\tanh(x)$')

    plt.xticks((-np.pi, 0.0, np.pi), (r'$-\pi$', '0', r'$\pi$'))
    plt.yticks((-1.0, 0.0, 1.0))

    plt.hlines(0.0, np.min(xx), np.max(xx), color='gray', linestyle='dashed')
    plt.plot(xx, tanh)
    plt.tight_layout()

    plt.savefig('tanh.svg', dpi=900, format='svg')
    plt.show()


def plot_sigmoid():
    '''Plots the sigmoid function in given interval.'''

    xx = np.linspace(-6.0, 6.0, 1000)
    sigmoid = 1.0 / (1.0 + np.exp(- xx))

    plt.figure(1, figsize=[7, 5])
    plt.title('Sigmoid Function')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$S(x)$')

    plt.yticks((0.0, 1.0))

    plt.hlines(0.5, np.min(xx), np.max(xx), color='gray', linestyle='dashed')
    plt.plot(xx, sigmoid)
    plt.tight_layout()

    plt.savefig('sigmoid.svg', dpi=900, format='svg')
    plt.show()


def plot_gaussian():
    '''Plots gaussian function in given interval.'''

    xx = np.linspace(-3.0, 3.0, 1000)
    gaussian = np.exp(- xx**2)

    plt.figure(1, figsize=[7, 5])
    plt.title('Gaussian Function')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$S(x)$')

    plt.yticks((0.0, 1.0))

    plt.vlines(0.0, 0.0, 1.0, color='gray', linestyle='dashed')
    plt.plot(xx, gaussian)
    plt.tight_layout()

    plt.savefig('gaussian.svg', dpi=900, format='svg')
    plt.show()


def plot_linear():
    '''Plots linear function in given interval.'''

    xx = np.linspace(-1.0, 1.0, 1000)
    linear = xx

    plt.figure(1, figsize=[7, 5])
    plt.title('Linear Function')
    plt.xlabel(r'$y$')
    plt.ylabel(r'$x$')

    plt.xticks((-1.0, 0.0, 1.0))
    plt.yticks((-1.0, 0.0, 1.0))

    plt.hlines(0.0, np.min(xx), np.max(xx), color='gray', linestyle='dashed')
    plt.plot(xx, linear)
    plt.tight_layout()

    plt.savefig('linear.svg', dpi=900, format='svg')
    plt.show()


def plot_relu():
    '''Plots relu function in given interval.'''

    xx = np.linspace(-1.0, 1.0, 1000)
    relu = xx * (xx > 0.0)

    plt.figure(1, figsize=[7, 5])
    plt.title('ReLU Function')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$R(x)$')

    plt.xticks((-1.0, 0.0, 1.0))
    plt.yticks((-1.0, 0.0, 1.0))

    plt.plot(xx, relu)
    plt.tight_layout()

    plt.savefig('relu.svg', dpi=900, format='svg')
    plt.show()


def plot_heaviside():
    '''Plots heaviside function in given interval.'''

    xx = np.linspace(-1.0, 1.0, 1000)
    heaviside = np.heaviside(xx, 0.5)

    plt.figure(1, figsize=[7, 5])
    plt.title('Heaviside Function')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$\Theta(x)$')

    plt.xticks((-1.0, 0.0, 1.0))
    plt.yticks((0.0, 0.5, 1.0))

    plt.plot(xx, heaviside)
    plt.tight_layout()

    plt.savefig('heaviside.svg', dpi=900, format='svg')
    plt.show()


if __name__ == '__main__':
    main()
