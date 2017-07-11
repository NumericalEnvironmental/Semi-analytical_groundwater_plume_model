# Semi-analytical_groundwater_plume_model
A julia-based script to compute groundwater plumes by superposition and numerical integration

This julia script integrates groundwater plume source terms to develop an advective-dispersive model for informing allocation calculations. Multiple sources with different shapes and time-dependent release histories are supported, as are scale-dependent dispersion coefficients. The metholodology is discussed in more detail on my blog (link pending). 

The following text input files are required:

* aquifer.txt - aquifer properties (groundwater velocity, flow direction, porosity, thickness, ratio of dispersivities to plume length, retardation coefficient, decay coefficient)
* grid.txt - grid maxima and minima, and discretization (for generating contour plots with external packages)
* sources.txt - source terms, including reference locations (in absolute coordinates), polynomials that describe mass influx as a function of time, source area start and stop extents (relative coordinates, with respect to groundwater flow direction), and polynomials that describe source area extents perpendicular to groundwater flow direction

An equivalent python script, which uses the same input files, is posted directly in the blog. The python script requires the numba just-in-time compiler from Continuum Analytics for improved performance.

I'd appreciate hearing back from you if you find the code useful. Questions or comments are welcome at walt.mcnab@gmail.com.

THIS CODE/SOFTWARE IS PROVIDED IN SOURCE OR BINARY FORM "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
