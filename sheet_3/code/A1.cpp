#include <eigen3/Eigen/Dense>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <vector>
#include <tuple>

using namespace std;
using namespace Eigen;

// =================================================================================================
//                      PROGRAM STRUCTURE
//
//                          ===========
//                          | Dataset |
//                          ===========
//                              |
//                           ========
//                           | Data |
//                           ========
//                              |
// ==============             ======                =============
// | Thermostat | ----------- | MD | -------------- | Potential |
// ==============             ======                =============
//                              |
//                           ========
//                           | main |
//                           ========
//
// - The MD class contains the primary logic of the algorithm
// - Thermostat and Potential class are separated from it to allow different
// thermostats and potentials can be implemented by inheritance and used flexibly.
// can be used
// - The Data class stores the data stored in the MD simulation and
// takes care of the storage
// - Dataset is a data set consisting of time, temperature, ...
// - Data holds several datasets and some more data that are not time-resolved
// stored
// - Instead of using getter and setter the members of Data
// and Dataset are public, since they are simple data containers.
// - main() calls MD with the parameters needed for the task parts
//
// Notes on the vectors used:
// - For performance reasons, we use Vector2d instead of VectorXd.
// However, this makes the possible transition to 3d states more cumbersome than with VectorXd
// - For lists of data std::vector is used.
// =================================================================================================

// CONSTANTS
double epsilon = 1;
double sigma = 1;

// ================================ Potential-class ================================================

// Virtual class from which concrete potentials can be inherited
// (Here only Lennard-Jones necessary, but so you can quickly implement other potentials).
class Potential {
  public:
    virtual double V(double r2) const = 0;    // Virtual function
    virtual Vector2d F(Vector2d r) const = 0; // Virtual function
};

class PotentialLJ : public Potential {
    // LJ = Lennard-Jones
  public:
    double V(double r2) const;    // Overwrites virtual function
    Vector2d F(Vector2d r) const; // Overwrites virtual function
};

double PotentialLJ::V(double r2) const {
    // For the potential, the square of the vector length is sufficient, which saves a root calculation.
    double r2_6 = pow(sigma / r2, 6);
    // cout << r2 << ", ";
    return 4 * epsilon * (r2_6 * r2_6 - r2_6);
}

Vector2d PotentialLJ::F(Vector2d r) const {
    // → siehe Kierfeld, S. 72

    double K = 24 * (-pow(r.norm(), -7.) + 2 * pow(r.norm(), -13.));
    return r.normalized() * K;
}

// ------------------------------ End of Potential-class -------------------------------------------
// ================================ Thermostat class ===============================================

// Virtual class from which concrete thermostats can be inherited
class Thermostat {
  public:
    virtual void rescale(vector<Vector2d> &v, double T, double T0) const = 0;
};

// No thermostat
class NoThermostat : public Thermostat {
  public:
    void rescale(vector<Vector2d> &v, double T, double T0) const {} // does nothing
};

// Isokinetic thermostat for task d)
class IsokinThermostat : public Thermostat {
  public:
    void rescale(vector<Vector2d> &v, double T, double T0) const;
};

void IsokinThermostat::rescale(vector<Vector2d> &v, double T, double T0) const {
    // Scale velocities to achieve desired temperature
    // NOTE: T ~ Ekin ~ p^2 ~ v^2 → v ~ sqrt(T)
    double scale = sqrt(T0 / T);
    for (int i = 0; i < v.size(); i++) {
        v[i] *= scale;
    }
    cout << "velocities scaled" << endl;
}

// ------------------------------ End of Thermostat class ------------------------------------------
// ================================ Data-Structs ===============================================

// Data set for time resolved data
// (structs basically the same as class, but all members are public by default)
struct Dataset {
    double t, T, Ekin, Epot;
    Vector2d vS;
};

// Return data of the MD simulation
// Data data(n) constructor; reserves memory and fills pair correlation function with 0s
struct Data {
    vector<Dataset> datasets; // Time-resolved datasets.
    vector<double> rBin, g;   // Averaged pair correlation function: g(rBin) (?)
    vector<Vector2d> r;       // snapshot of the final position
                              // For task e) it may be useful to use r instead
                              // in the time-resolved datasets instead

    Data(uint n, uint numBins, double binSize);
    void save(const string &filenameSets,
              const string &filenameG,
              const string &filenameR) const;
};

Data::Data(uint n, uint numBins, double binSize)
    : datasets(n), // Initializer list, because it calls constructors of the members
      rBin(numBins),
      g(numBins, 0.),
      r(0) {
}

void Data::save(const string &filenameSets, const string &filenameG, const string &filenameR) const {
    ofstream myfile;

    // save datasets
    myfile.open(filenameSets);
    myfile << "t\tT\tEkin\tEpot\tvSx\tvSy" << endl;

    for (const Dataset &set : datasets) {
        myfile << set.t << "\t" << set.T << "\t" << set.Ekin << "\t" << set.Epot << "\t" << set.vS.x() << "\t" << set.vS.y() << endl;
    }
    myfile.close();

    // save (final) pair correlation function
    myfile.open(filenameG);
    myfile << "rBin\tg" << endl;
    for (int i = 0; i < g.size(); i++) {
        myfile << rBin[i] << "\t" << g[i] << endl;
    }
    myfile.close();

    // save final positions
    myfile.open(filenameR);
    myfile << "x\ty" << endl;
    for (int i = 0; i < r.size(); i++) {
        myfile << r[i].x() << "\t" << r[i].y() << endl;
    }
    myfile.close();
}

// ------------------------------ End of Data-Structs ------------------------------------------
// ================================ MD-Class ===============================================

class MD {
  public:
    MD(double L, uint N, uint particlesPerRow, double T,
       Potential &potential, Thermostat &thermostat,
       uint numBins = 1000);

    void equilibrate(const double dt, const unsigned int n);
    Data measure(const double dt, const unsigned int n);

  private:
    void verlet(const double dt);
    tuple<vector<double>, vector<double>> calcPairCorrelation();

    vector<Vector2d> r, v;
    double L;
    uint N;
    Potential &potential;
    Thermostat &thermostat;
    double t = 0.;

    uint numBins;
    double binSize;

    double T0; // ours

    // Particles are moved in box [0,L]x[0,L].
    void centerParticles();

    // Calculations of important measured variables
    double calcT() const;
    double calcEkin() const;
    double calcEpot() const;
    Vector2d calcvS() const;
    Dataset calcDataset() const;

    // Calculation of the acceleration
    // To avoid redundant calculations, it may be useful to update the histogram
    // when calculating the accelerations, so it is passed here as a reference
    vector<Vector2d> calcAcc(vector<double> &hist) const; // not used
    vector<Vector2d> calcAcc_2(vector<double> &hist) const;

    // Calculation of the distance vector between particle r[i] and closest mirror particle of r[j]
    Vector2d calcDistanceVec(uint i, uint j) const;
};

// Initialization of the system via constructor
MD::MD(double L, uint N, uint particlesPerRow, double T0,
       Potential &potential, Thermostat &thermostat,
       uint numBins)
    : L(L),
      N(N),
      potential(potential),
      thermostat(thermostat),
      numBins(numBins),
      // /2 = berücksichtige cutoff…
      binSize(L / numBins / 2), // TODO
      T0(T0) {

    cout << "Bin size: " << binSize << endl;
    cout << "L: " << L << endl;
    cout << "numBins: " << numBins << endl;

    mt19937 rnd;
    uniform_real_distribution<double> dist(0, 1);

    // NOTE: distance between particles: L / particlesPerRow = 2*n*sigma / n = 2*sigma
    for (uint i = 0; i < particlesPerRow; i++) {
        for (uint j = 0; j < particlesPerRow; j++) {
            // NOTE: We offset by 1 so the particles are not on the boundary
            Vector2d r_ij = Vector2d(2 * i * sigma + 1, 2 * j * sigma + 1);
            r.push_back(r_ij);

            Vector2d v_ij = Vector2d(dist(rnd), dist(rnd));
            v.push_back(v_ij);
        }
    }
    cout << "particles created" << endl;

    // subtract mean velocity (vS) from all velocities
    Vector2d vS = calcvS();
    for (int i = 0; i < N; i++) {
        v[i] -= vS;
    }
    cout << "mean velocity subtracted" << endl;

    // Scale velocities to achieve desired temperature
    // NOTE: T ~ Ekin ~ p^2 ~ v^2 → v ~ sqrt(T)
    double T = calcT();
    double scale = sqrt(T0 / T);
    for (int i = 0; i < N; i++) {
        v[i] *= scale;
    }
    cout << "velocities scaled" << endl;

}

// Integration without data acquisition for pure equilibration
void MD::equilibrate(const double dt, const unsigned int n) {
    vector<double> stupidHist;
    cout << "Equilibrating...";
    for (uint i = 1; i <= n; i++) {
        verlet(dt);
        cout << ".";
    }
    cout << endl;
}

Data MD::measure(const double dt, const unsigned int n) { // n=steps???
    // double t, T, Ekin, Epot;
    // Vector2d vS;

    // create data object and save initial values
    Data data(n + 1, numBins, binSize); // NOTE: n+1 because we want to save the initial values too
    data.datasets[0] = calcDataset();
    data.r = r;
    vector<vector<double>> g_set(n, vector<double>(numBins, 0.)); // save g for each timestep

    ofstream file("build/r_merged.txt");

    // verlet algorithm
    for (uint i = 1; i <= n; i++) {
        t = i * dt;
        verlet(dt);

        // write r to csv file
        for (uint j = 0; j < N; j++) {
            file << r[j].x() << "," << r[j].y() << " | ";
        }
        file << endl;

        // save data in dataset
        data.datasets[i] = calcDataset();
        data.r = r;

        // update pair correlation function
        // form the mean value of all pair correlation functions

        tie(data.rBin, g_set[i-1]) = calcPairCorrelation();
    }

    // time average over g_set
    for (int i = 0; i < n; i++) {
        for (int binIndex = 0; binIndex < numBins; binIndex++) {
            data.g[binIndex] += g_set[i][binIndex] / n;
        }
    }

    file.close();

    return data;
}

void MD::verlet(const double dt) {
    vector<double> stupidHist; // dummy

    vector<Vector2d> a_n = calcAcc_2(stupidHist); // calc acceleration for each particle
    vector<Vector2d> r_n = r;                     // save r_n

    // calc r_n+1 for each particle
    for (uint ind = 0; ind < N; ind++) {
        // r_n+1 = r_n + v_n*dt + 0.5*a_n*dt*dt
        r[ind] = r_n[ind] + v[ind] * dt + 0.5 * a_n[ind] * dt * dt;
    }

    // periodic boundary conditions
    centerParticles();

    // calc acceleration for new positions: a_n+1
    vector<Vector2d> a_np1 = calcAcc_2(stupidHist);

    // calc v_n+1 for each particle
    for (uint ind = 0; ind < N; ind++) {
        // v_n+1 = v_n + 0.5*(a_n+1 + a_n)*dt
        v[ind] = v[ind] + 0.5 * (a_np1[ind] + a_n[ind])*dt;
    }

    // rescale velocities
    thermostat.rescale(v, calcT(), T0); // TODO: reuse T from calcDataset()?
}

// Particles are moved in box [0,L]x[0,L].
void MD::centerParticles() {
    // Ensure periodic boundary conditions ([0,L]x[0,L])
    // If the particle is outside of box, move it to the other side
    for (int i = 0; i < N; i++) {
        r[i].x() = fmod(fmod(r[i].x(), L) + L, L);
        r[i].y() = fmod(fmod(r[i].y(), L) + L, L);
    }
}

double MD::calcT() const {
    double N_f = 3 * N - 3; // number of degrees of freedom

    return 2 * calcEkin() / N_f; // TODO: epsilon, k_B?
}

double MD::calcEkin() const {
    double Ekin = 0;
    for (int i = 0; i < N; i++) {
        Ekin += 0.5 * v[i].squaredNorm();
    }
    return Ekin;
}

double MD::calcEpot() const {
    double Epot = 0;
    // cout << "BEGIN V_ij: ";
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            Vector2d r_ij = calcDistanceVec(i, j);
            double V_ij = potential.V(r_ij.squaredNorm());

            // if r_ij is (0,0), then the particles are too far away from each other
            if (r_ij.squaredNorm() == 0) {
                continue;
            }

            // cout << V_ij << ", ";
            Epot += V_ij;
        }
    }
    return Epot;
}

Vector2d MD::calcvS() const {
    // vS means Schwerpunktsgeschwindigkeit
    Vector2d vS(0, 0);
    for (int i = 0; i < N; i++) {
        // NOTE: equal mass particles
        vS += v[i];
    }
    return vS / N;
}

Dataset MD::calcDataset() const {
    // calc data
    double T = calcT();       // calc temperature
    double Ekin = calcEkin(); // calc energy
    double Epot = calcEpot(); // calc Epot
    Vector2d vS = calcvS();   // calc vS

    Dataset set = {t, T, Ekin, Epot, vS};

    cout << "\nTime: " << t << endl;
    cout << "T: " << T << ", Ekin: " << Ekin << ", Epot: " << Epot << ", vS: (" << vS.x() << ", " << vS.y() << ")" << endl;

    return set;
}

Vector2d MD::calcDistanceVec(uint i, uint j) const {
    // return r[i] - r[j];

    const double r_cutoff = L / 2;

    // boundary conditions and cutoff
    for (int x = -1; x <= 1; x++) {
        for (int y = -1; y <= 1; y++) {
            Vector2d r_ij = (r[i] - r[j]) + Vector2d(x * L, y * L);
            if (r_ij.squaredNorm() < r_cutoff * r_cutoff) {
                return r_ij;
            }
        }
    }

    return Vector2d(0, 0); // FIXME: suboptimal
}

vector<Vector2d> MD::calcAcc(vector<double> &hist) const {
    vector<Vector2d> a(N, Vector2d(0, 0));

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < i; j++) { // NOTE: j < i, because F_ij = -F_ji
            // cout << "i: " << i << ", j: " << j << endl;
            Vector2d r_ij = calcDistanceVec(i, j);

            // if r_ij is (0,0), then the particles are too far away from each other
            if (r_ij.squaredNorm() == 0) {
                continue;
            }

            Vector2d F_ij = potential.F(r_ij); // NOTE: F = m*a and m = 1
            a[i] += F_ij;
            a[j] -= F_ij;
        }
    }

    return a;
}

vector<Vector2d> MD::calcAcc_2(vector<double> &hist) const {
    // variant of calcAcc that uses explicit for loops
    vector<Vector2d> a(N, Vector2d(0, 0));

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            // self-interaction is not allowed
            // NOTE: This would be caught by the if statement below, but this is more efficient
            if (i == j) {
                continue;
            }

            Vector2d r_ij = calcDistanceVec(i, j);
            Vector2d F_ij = potential.F(r_ij);

            // if r_ij is (0,0), then the particles are too far away from each other
            if (r_ij.squaredNorm() == 0) {
                continue;
            }

            a[i] += F_ij;
        }
    }
    return a;
}

tuple<vector<double>, vector<double>> MD::calcPairCorrelation() {
    // vector<double> rBin, g
    // numBins, binSize

    // iterate over particles
    //  → iterate over other particles
    //     → calc distance
    //     → calc bin index
    //     → add to g
    // average over these pair correlation functions

    vector<double> rBin_temp(numBins, 0.);
    for (int i = 1; i <= numBins; i++) {
        rBin_temp[i] = i * binSize;
    }
    // NOTE: r will denote the right edge of the bin! ⚠️

    // initialize temporary g with "shape" (N, numBins)
    vector<vector<double>> g_i(N, vector<double>(numBins, 0.));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double r_ij = calcDistanceVec(i, j).norm();

            if (i == j || r_ij == 0) {
                // skip if particle is itself or beyond cutoff
                continue;
            }

            int binIndex = r_ij / binSize;
            g_i[i][binIndex] += 1;
        }
    }

    vector<double> g_temp(numBins, 0.);
    // average over g_i
    for (int i = 0; i < N; i++) {
        for (int binIndex = 0; binIndex < numBins; binIndex++) {
            g_temp[binIndex] += g_i[i][binIndex] / N;
        }
    }

    // weight by bin size (area)
    for (int binIndex = 0; binIndex < numBins; binIndex++) {
        double r = binIndex * binSize;
        double A = M_PI * (pow(r + binSize, 2) - pow(r, 2));
        g_temp[binIndex] /= A;
    }

    // norm to 1
    double norm = 0;
    for (int binIndex = 0; binIndex < numBins; binIndex++) {
        norm += g_temp[binIndex];
    }
    for (int binIndex = 0; binIndex < numBins; binIndex++) {
        g_temp[binIndex] /= norm;
    }

    // return a tuple of r and g
    return make_tuple(rBin_temp, g_temp);
}

// ------------------------------ End of MD-class ------------------------------------------

int main(void) {
    cout << "main start" << endl;

    PotentialLJ LJ;
    NoThermostat noThermo;
    IsokinThermostat isoThermo;

    // PARAMETERS
    const uint particlesPerRow = 8;
    const uint N = particlesPerRow * particlesPerRow;
    const double L = 2 * particlesPerRow * sigma;
    const int numBins = 500;

    // b) Equilibration test (example for T=1, dt=0.01, steps_meas=1000)
    double T_b = 1;     // T(0)
    double dt_b = 0.01;
    uint steps_b = 1000;

    MD md_b(L, N, particlesPerRow, T_b, LJ, noThermo, numBins);
    cout << "++ MD init complete" << endl;

    md_b.measure(dt_b, steps_b).save("build/b)set.tsv", "build/b)g.tsv", "build/b)r.tsv");
    cout << "++ MD measure complete" << endl;

    // c) Pair correlation function
    string TstringVec[3] = {"0.01", "1", "100"};
    string dtstringVec[3] = {"0.01", "0.01", "0.001"};

    uint equiSteps_c = 2000;
    uint steps_c = 10000;

    for (int i = 0; i < 3; i++) {
        double T_c = stod(TstringVec[i]);
        double dt_c = stod(dtstringVec[i]);
        
        MD md_c(L, N, particlesPerRow, T_c, LJ, noThermo, numBins);
        cout << "++ MD equilibrate start for T = " << T_c << " and dt = " << dt_c << endl;
        md_c.equilibrate(dt_c, equiSteps_c);

        cout << "++ MD measure start for T = " << T_c << " and dt = " << dt_c << endl;
        md_c.measure(dt_c, steps_c).save("build/c)set_" + TstringVec[i] + ".tsv", "build/c)g_" + TstringVec[i] + ".tsv", "build/c)r_" + TstringVec[i] + ".tsv");
    }

    // d) Isokinetic thermostat
    double T_d = 0.01;
    double dt_d = 0.01;       
    uint equiSteps_d = 5000;
    uint steps_d = 10000;

    MD md_d = MD(L, N, particlesPerRow, T_d, LJ, isoThermo, numBins);

    md_d.equilibrate(dt_d, equiSteps_d);
    md_d.measure(dt_d, steps_d).save("build/d)set.tsv", "build/d)g.tsv", "build/d)r.tsv");

    return 0;
}