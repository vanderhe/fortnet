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
    plot_atan()
    plot_sigmoid()
    plot_softplus()
    plot_gaussian()
    plot_linear()
    plot_relu()
    plot_lrelu()
    plot_bent()
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


def plot_atan():
    '''Plots the arcus tangent in given interval.'''

    xx = np.linspace(-10.0, 10.0, 1000)
    atan = np.arctan(xx)

    plt.figure(1, figsize=[7, 5])
    plt.title('Arcus Tangent')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$\arctan(x)$')

    plt.xticks((-10.0, 0.0, 10.0))
    plt.ylim((-np.pi / 2.0, np.pi / 2.0))
    plt.yticks((-np.pi / 2.0, 0.0, np.pi / 2.0),
               (r'$-\frac{\pi}{2}$', '0', r'$\frac{\pi}{2}$'))

    plt.hlines(0.0, np.min(xx), np.max(xx), color='gray', linestyle='dashed')
    plt.plot(xx, atan)
    plt.tight_layout()

    plt.savefig('atan.svg', dpi=900, format='svg')
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


def plot_softplus():
    '''Plots the SoftPlus function in given interval.'''

    xx = np.linspace(-3.0, 3.0, 1000)
    softplus = np.log(1.0 + np.exp(xx))

    plt.figure(1, figsize=[7, 5])
    plt.title('SoftPlus Function')
    plt.xlabel(r'$x$')
    plt.ylabel(r'SoftPlus')

    plt.yticks((-1.0, 0.0, 1.0, 2.0, 3.0))

    plt.plot(xx, softplus)
    plt.tight_layout()

    plt.savefig('softplus.svg', dpi=900, format='svg')
    plt.show()


def plot_gaussian():
    '''Plots gaussian function in given interval.'''

    xx = np.linspace(-3.0, 3.0, 1000)
    gaussian = np.exp(- xx**2)

    plt.figure(1, figsize=[7, 5])
    plt.title('Gaussian Function')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$G(x)$')

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
    '''Plots ReLU function in given interval.'''

    xx = np.linspace(-1.0, 1.0, 1000)
    relu = xx * (xx > 0.0)

    plt.figure(1, figsize=[7, 5])
    plt.title('ReLU Function')
    plt.xlabel(r'$x$')
    plt.ylabel('ReLU')

    plt.xticks((-1.0, 0.0, 1.0))
    plt.yticks((-1.0, 0.0, 1.0))

    plt.plot(xx, relu)
    plt.tight_layout()

    plt.savefig('relu.svg', dpi=900, format='svg')
    plt.show()


def plot_lrelu():
    '''Plots leaky ReLU function in given interval.'''

    xx = np.linspace(-1.0, 1.0, 1000)
    lrelu = np.array([max(0.01 * val, val) for val in xx], dtype=float)

    plt.figure(1, figsize=[7, 5])
    plt.title('Leaky ReLU Function')
    plt.xlabel(r'$x$')
    plt.ylabel('Leaky ReLU')

    plt.xticks((-1.0, 0.0, 1.0))
    plt.yticks((-1.0, 0.0, 1.0))

    plt.hlines(0.0, np.min(xx), np.max(xx), color='gray', linestyle='dashed')

    plt.plot(xx, lrelu)
    plt.tight_layout()

    plt.savefig('lrelu.svg', dpi=900, format='svg')
    plt.show()


def plot_bent():
    '''Plots Bent identity function in given interval.'''

    xx = np.linspace(-1.0, 1.0, 1000)
    bent = (np.sqrt(xx**2 + 1.0) - 1.0) / 2.0 + xx

    plt.figure(1, figsize=[7, 5])
    plt.title('Bent Identity Function')
    plt.xlabel(r'$x$')
    plt.ylabel('Bent identity')

    plt.xticks((-2.0, -1.0, 0.0, 1.0))
    plt.yticks((-1.0, 0.0, 1.0))

    plt.hlines(0.0, np.min(xx), np.max(xx), color='gray', linestyle='dashed')

    plt.plot(xx, bent)
    plt.tight_layout()

    plt.savefig('bent.svg', dpi=900, format='svg')
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
