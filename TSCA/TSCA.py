from warnings import warn
from scipy import *
from scipy import signal
from pylab import *

def make_autocorrelation_matrix(signal1):
    mat = signal1[:,newaxis] * signal1[newaxis,:]
    mat /= mat.max()
    return mat

# TODO read up on Weiner-Knichne . . . 
class TSCA:
    def __init__(self, image_temporal_data):
        assert image_temporal_data.ndim == 3, "Data must be three-dimensional (x by y by t)"
        self.shape = image_temporal_data.shape[0], image_temporal_data.shape[1]
        temporal_data = image_temporal_data.reshape((self.shape[0] * self.shape[1], -1))
        self.Z = temporal_data

        self.channel_count, self.time_points = self.Z.shape

        self.temporal_signals = []
        self.temporal_matrices = []
        self.signal_count = 0

        # Status Variables
        self.mat_UTD = False

    def add_signal(self, signal_structure):
        assert signal_structure.shape[0] == self.time_points, "Time signal must match time course length"
        self.temporal_signals  += [ signal_structure ]
        self.temporal_matrices += [ make_autocorrelation_matrix(signal_structure) ]
        self.signal_count += 1
        self.mat_UTD = False
        
    def remove_signal(self, number):
        assert number < self.signal_count, "Invalid signal number . . ."
        # TODO, optimize (throw out matrices instead of recomputing all)
        del self.temporal_matrices[number]
        self.signal_count -= 1
        self.mat_UTD = False

    def compute_matrices(self, gammae):
        assert self.signal_count > 0, "Add signals first!"
        assert gammae.shape[0]==len(self.temporal_signals), "You mush have one weight per signal"

        print "\t Computing cross correlation matrices"
        # Compute the Correlation and Cross-Correlation terms
        self.intermediate_matrix = matrix( zeros((self.signal_count, self.signal_count)) )
        for n1, m1 in enumerate(self.temporal_matrices):
            for n2, m2 in enumerate(self.temporal_matrices):
                self.intermediate_matrix[n1, n2] = sum(m1 * m2)
        self.intermediate_matrix = self.intermediate_matrix.I

        # Weightings and RHS of equation TODO (number from paper)
        traces = []
        for n1, m1 in enumerate(self.temporal_matrices):
            traces += [trace(m1)]
        
        # Compute the traces to norm the covariance projections
        traces = array(traces)
        self.weights = traces * gammae
        self.q_mat_comp = array(dot(self.intermediate_matrix, self.weights)).flatten()

        print "Computing Q"
        # Make the matrix Q described by eq. TODO (number)
        self.Q = matrix( zeros((self.time_points, self.time_points)) )
        for comp, mat  in enumerate(self.temporal_matrices):
            print mat.shape,  self.q_mat_comp[comp].shape
            i = self.q_mat_comp[comp] * mat
            self.Q += i
        
        # Create matrix that we solve to find phi
        self.ZQZT = self.Z * self.Q * self.Z.T

        # Update "up to date" flag
        self.mat_UTD = True
        
    def solve(self):
        # If the matrices are configed, weightings have alredy been specified
        if not self.mat_UTD:
            warn("Weightings of components not specified.  Using all equal!")
            self.compute_matrices( ones(self.signal_count) )

        print "Computing Eigenvectors of", self.ZQZT.shape, "matrix . . ."
        # self.evals, self.evec = eigh(self.ZQZT)
        self.evals, self.evec = eig(self.ZQZT)
        self.valorder = self.evals.argsort()

    def eigenimage(self, phi_num):
        phi = self.evec[:, self.valorder[phi_num]]
        phi = real( phi.reshape(self.shape) )
        return phi


if __name__ == "__main__":
    print "Running Test Case Similar to TSCA paper!"
    size_of_image = (30, 30)
    time_slices = arange(0,1000)

    x, y = mgrid[-1.0:1.0:size_of_image[0]*1j, -1.0:1.0:size_of_image[1]*1j]

    # The Geometric Shapes
    in_circle = 1.0 * (sqrt(x**2 + y**2) < 0.5)
    x_grating = cos(x * 2 * pi)

    # Their individual time Courses
    circ_strength = cos(time_slices * 0.21342 / (2 * pi))
    wave_strength = cos(time_slices * 0.13420 / (2 * pi))

    # X x Y x T
    time_course = in_circle[:,:,newaxis] * circ_strength + x_grating[:,:,newaxis] * wave_strength

    ts = TSCA(time_course)
    ts.add_signal(circ_strength)
    ts.add_signal(wave_strength)
    ts.solve()
