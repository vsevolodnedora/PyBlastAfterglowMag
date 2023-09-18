//
// Created by vsevolod on 21/12/22.
//

#ifndef SRC_UTILS_H
#define SRC_UTILS_H

#include "pch.h"
#include "logger.h"

//typedef std::valarray<double> Array;
//typedef std::vector<Array> VecArray;
typedef std::vector<double> Vector;
typedef std::vector<std::vector<double>> VecVector;
typedef std::unordered_map<std::string,double> StrDbMap;
typedef std::unordered_map<std::string,std::string> StrStrMap;

inline namespace CGS{
    const double &c  = 2.99792458e10;     // speed of light in cm / s
    const double &pi = 3.141592653589793; // pi
    const double &mp = 1.6726e-24;        // proton mass
    const double &me = 9.1094e-28;        //electron mass
    const double &mppme = mp + me;
    const double &qe = 4.803204e-10;
    const double &sigmaT = 6.6524e-25;
    const double &gamma_c_w_fac = 6 * pi * me * c / sigmaT; // Used in gamma_c calculation
    const double &kB = 1.380658e-16; // Boltzmann constant
    const double &hcgs = 6.6260755e-27; // Planck constant in cgs
    const double &h    = 6.6260755e-27; // erg s
    const double &lambda_c = (hcgs / (me * c)); // Compton wavelength
    const double &mec2 = 8.187105649650028e-07;  // erg # electron mass in erg, mass_energy equivalence
    const double &pc = 3.0857e18; // cm
    const double &year= 3.154e+7; // sec
    const double &day = 86400;
    const double &solar_m = 1.989e+33; // solar mass in g
    const double &deg2rad = 0.017453292;
    const double &cgs2mJy = 1.e26;
    const double &gravconst = 6.67259e-8;// # cm^3 g^-1 s^-2
    /// --------------------------------
    const double &SI_pc   =3.0856E+16; // m
    const double &RtD = 57.295779513082325;
    const double &A_RAD = 7.565767e-15; // radiation constant in unit of [erg/cm^3/K^4]
    const double &EV_TO_ERG = 1.602e-12; // [erg/eV]
    const double &M_ELEC = 9.11e-28; // electron mass in unit of [g]
    const double &ELEC = 4.803e-10; // electron charge in unit of [cgs]
    const double &SIGMA_T = 6.65e-25; // Thomson cross section in unit of [cm^2]
    const double &K_B = 1.38e-16; // Boltzman constant in unit of [erg/K]
    const double &H = 6.626e-27; // Planck constant in unit of [erg s]
    const double &M_PRO = 1.0/6.02e23; /* atomic mass in unit of [g] */
    const double &MeC2 = 8.1871398e-7; // in unit of [erg]
    const double &GAMMA13 = 2.67893; /* Gamma(1/3) */
    const double &rad2mas = 206264806.247; // radians to milli-arcseconds
};

namespace TOOLS{
    // --- make a vector of logarithmically spaced points
    template<typename T>
    class Logspace {
    private:
        T curValue, base;

    public:
        Logspace(T first, T base) : curValue(first), base(base) {}

        T operator()() {
            T retval = curValue;
            curValue *= base;
            return retval;
        }
    };

    static void MakeLogspaceVec(std::vector<double> & vec, const double &start, const double &stop, const double &base = 10) {
        int num = (int)vec.size();
        double realStart = std::pow(base, start);
        double realBase = std::pow(base, (stop-start)/num);
        std::generate_n(std::back_inserter(vec), num, Logspace<double>(realStart,realBase));
    }

    static std::vector<double> MakeLogspaceVec(const double &start,
                                               const double &stop,
                                               const int &num = 50,
                                               const double &base = 10) {
        double realStart = pow(base, start);
        double realBase = pow(base, (stop-start)/num);

        std::vector<double> retval;
        retval.reserve(num);

        std::generate_n(std::back_inserter(retval), num, Logspace<double>(realStart,realBase));
        return std::move(retval);

    }

    //
//    static std::valarray<double> MakeLogspace(const double &start,
//                                              const double &stop,
//                                              const int &num = 50,
//                                              const double &base = 10) {
//        auto retval = MakeLogspaceVec(start, stop, num, base);
//        return std::move(std::valarray<double> (retval.data(), retval.size()));
//
//    }

    std::vector<double> linspace(double first, double last, int len) {
        std::vector<double> result(len);
        double step = (last-first) / (len - 1);
        for (int i=0; i<len; i++) { result[i] = first + i*step; }
        return result;
    }

}

struct Timer {

    //std::chrono::time_point<std::chrono::steady_clock> start, end;
    std::chrono::system_clock::time_point start, end;
    std::chrono::duration<float> duration;

    Timer() {
        start = std::chrono::high_resolution_clock::now();
    }

    auto checkPoint() const{
        using namespace std::literals::chrono_literals;
        auto m_end = std::chrono::high_resolution_clock::now();
        auto m_duration = m_end-start; // duration in seconds
        auto duration_ms = m_duration.count()/1.e9;// * 1000.0f;
        return duration_ms; // s ???
//        std::cout << "Time took " << duration_ms << " ms" << std::endl;
    }

//    ~Timer()
//    {
//        end = std::chrono::high_resolution_clock::now();
//        duration = end-start; // duration in seconds
//        float duration_ms = duration.count() * 1000.0f;
//        std::cout << "Time took " << duration_ms << " ms" << std::endl;
//    }

};

double maxValue( const Vector & vec )
{
    double val_max = vec[0]; // What could go wrong on this line?

    for ( int i = 0; i < vec.size(); ++i ) // Is the first element necessary?
    {
        if ( val_max < vec[i] )
        {
            val_max = vec[i];
        }
    }
    return val_max;
}

double findMinimum(Vector & vec) {
    double minimum = std::numeric_limits<double>::max();

    for (const double& value : vec) {
        if (value < minimum) {
            minimum = value;
        }
    }
    return minimum;
}
double findMaximum(const Vector & vec) {
    double maximum = std::numeric_limits<double>::min();

    for (const double& value : vec) {
        if (value > maximum) {
            maximum = value;
        }
    }

    return maximum;
}

/// https://stackoverflow.com/questions/42533070/finding-minimum-element-of-a-vector-in-c
size_t indexOfMinimumElement(const std::vector<double>& input)
{
    if (input.empty())
        return -1;
    auto ptrMinElement = std::min_element(input.begin(), input.end());
    return std::distance(input.begin(), ptrMinElement);
}

bool getBoolOpt(std::string opt, StrStrMap & opts,
                const char *loc, std::unique_ptr<logger> & p_log, const bool def, const bool req){
//    opt = "use_dens_prof_behind_jet_for_ejecta";
    bool b_val;
    std::string val;
    if ( opts.find(opt) == opts.end() )
        if (req){
//            std::cerr << AT << " \n Required parameter " << opt << " not found. exiting...\n";
            (*p_log)(LOG_ERR, loc) << "Required parameter " << opt << " not found. exiting...\n";
            exit(1);
        }
        else{
            (*p_log)(LOG_WARN, loc) << "Option " << opt << " not given. Using default value="<<def<<"\n";
            b_val = def;
        }
    else{
        val = (std::string)opts.at(opt);
        if (val == "yes"){
            b_val = true;
        }
        else if (val == "no"){
            b_val = false;
        }
        else{
            (*p_log)(LOG_ERR, loc) << " option for: " <<opt<<" given: " << opts.at(opt) << " is not recognized \n"
                                   << "Possible options: " << " yes " << " no " << "\n";
//            std::cerr << AT << " option for: " <<opt
//                      <<" given: " << opts.at(opt)
//                      << " is not recognized \n";
//            std::cerr << "Possible options: " << " yes " << " no " << "\n";
            exit(1);
        }
    }
    return b_val;
}
std::string getStrOpt(std::string opt, StrStrMap  opts,
                      const char *loc, std::unique_ptr<logger> & p_log, const std::string& def, const bool req){
//    opt = "use_dens_prof_behind_jet_for_ejecta";
//    std::string b_val;
    std::string val;
    if ( opts.find(opt) == opts.end() )
        if (req){
            (*p_log)(LOG_ERR, loc) << "Required parameter " << opt << " not found. exiting...\n";
//            std::cerr << AT << " \n Required parameter " << par << " not found. exiting...\n";
            exit(1);
        }
        else{
            (*p_log)(LOG_WARN, loc) << "Option " << opt << " not given. Using default value="<<def<<"\n";
            val = def;
        }
    else{
//        std::cout << "\t"<<val<<"\n";
        val = opts.at(opt);
    }


    return val;
//    throw std::runtime_error("error...");
}
double getDoublePar(std::string par, StrDbMap & pars,
                    const char *loc, std::unique_ptr<logger> & p_log, double def, const bool req){
    double val;
    if ( pars.find(par) == pars.end() ){
        if (req){
            (*p_log)(LOG_ERR, loc) << " Required parameter " << par << " not found. exiting...\n";
//            (*p_log)(LOG_ERR,loc) << " Required parameter " << par << " not found. exiting...\n";
            exit(1);
        }
        else{
            (*p_log)(LOG_WARN, loc) << "Parameter " << par << " not given. Using default value="<<def<<"\n";
            val = def;
        }
    }
    else{
        val = (double)pars.at(par);
    }
    return val;
}

Vector arrToVec(Vector & array){
    std::vector<double> vec;
    vec.assign(std::begin(array), std::end(array));
    return std::move( vec );
}
//void vecToArr(Vector & source, Array & target){
//    if ((target.size() != source.size())){
//        target.resizeEachImage(source.size(), 0.0);
//    }
//    for (size_t i = 0; i < source.size(); ++i)
//        target[i] = source[i];
//}


double linearExtrapolate(double x1, double x2, double y1, double y2, double new_x){
    double y;
    y = y1 + (new_x - x1) / (x2 - x1) * (y2 - y1);
    return y;
}

/// https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
template <typename T>
void sort_indexes(const std::vector<T> &v, std::vector<size_t> & idx) {

    // initialize original index locations
//    std::vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values
    stable_sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) { return v[i1] < v[i2]; });

//    return std::move( idx );
}
template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {

    // initialize original index locations
    std::vector<size_t> idx(v.size());
    std::iota(begin(idx), end(idx), 0);

    // sort indexes based on comparing values in vv
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when vv contains elements of equal values
    std::stable_sort(begin(idx), end(idx),
                     [&v](size_t i1, size_t i2)
                     {return v[i1] < v[i2];});

    return std::move( idx );
}
Vector sort_by_indexes(const Vector & array, const std::vector<size_t> & indexes){
    if (array.size() != indexes.size()){
        std::cerr << AT << " size mismatch\n";
        exit(1);
    }
    Vector sorted (array.size());
    for(size_t i = 0; i < array.size(); i++){
        sorted[i] = array[indexes[i]];
    }
    return std::move( sorted );
}

/**
 * Overload << for valarray and vector
 */
template<typename T>
std::ostream & operator<<(std::ostream & os, std::valarray<T> arr)
{
    os<<"{ ";
    std::copy(begin(arr), end(arr), std::ostream_iterator<T>(os, ", "));
    os<<"}";
    return os;
}
template<typename T>
std::ostream & operator<<(std::ostream & os, std::vector<T> arr)
{
    os<<"( ";
    std::copy(begin(arr), end(arr), std::ostream_iterator<T>(os, ", "));
    os<<")";
    return os;
}

/**
 * Format string aka printf()
 * @tparam Args
 * @param format
 * @param args
 * @return
 */
template<typename ... Args>
std::string string_format( const std::string& format, Args ... args )
{
    int size_s = std::snprintf( nullptr, 0, format.c_str(), args ... ) + 1; // Extra space for '\0'
    if( size_s <= 0 ){
        throw std::runtime_error( "Error during formatting." );
    }
    auto size = static_cast<size_t>( size_s );
    auto buf = std::make_unique<char[]>( size );
    std::snprintf( buf.get(), size, format.c_str(), args ... );
    return std::string( buf.get(), buf.get() + size - 1 ); // We don't want the '\0' inside
}

/*
 * Utils for saving data
 */
template<typename vecarr>
std::ostream& stream_arr(std::ostream& os, vecarr &data, std::vector<std::string> &names, std::unique_ptr<logger> & p_log)
{


    size_t n_arrs = 0;
    n_arrs = data.size();
    size_t n_vals = 0;
    n_vals = data[0].size();
    size_t n_names = 0;
    n_names = (const size_t)names.size();


    if (n_names != n_arrs){
        (*p_log)(LOG_ERR, AT) << "Error! Size of names_lc and of data vector is not the same "
                             "[n_names="<<n_names<<", n_arrs="<<n_arrs<<"]"
                          << "\n Exiting...\n";
//        std::cerr << AT<< "\n";
        exit(1);
    }
    if (data.isEmpty()) {
        (*p_log)(LOG_ERR,  AT) << "Error! Empty vector of data cannot be outputed"
                          << "Exiting...\n";
        std::cerr << AT<< "\n";
        exit(1);
    }
    if (data[0].size() == 0) {
        (*p_log)(LOG_ERR, AT) << "Error! Empty data array cannot be outputed"
                          << "\n Exiting...\n";
//        std::cerr << AT<< "\n";
        exit(1);
    }
    for (size_t j = 0 ; j < n_arrs; j ++){
        if (data[j].size() != data[0].size()){
            (*p_log)(LOG_WARN, AT) << "Error! Data arrays are of different length! " <<"\n"
                               << "       Size arr[0]=" << data[0].size() << " and arr[" << j << "]=" << data[j].size()
                               << "       For quantity j=" << j << " aka (" << names[j] << ")" << "\n";
//                      <<"Exiting...\n";
        }
    }

    std::vector<std::string> tmp;
    std::copy(names.begin(), names.end(), std::back_inserter(tmp));

    tmp.insert(tmp.begin(), "# ");

    // print the first row of var names_lc with '#' for (for PyBlastAfterglow read)
    for (size_t i = 0 ; i < n_names+1; i++){
        os << tmp[i]<<" ";
    }
    os<<"\n";

    // print the data, var by var
    for (size_t i = 0 ; i < n_vals; i++){
        for (size_t j = 0 ; j < n_arrs; j++){
            os << string_format("%.8e", data[j][i]) << " ";
        }
        os<<"\n";
    }
    //    for (auto & iarr : data) {
    //        for (int j = 0; j < iarr.size(); ++j) {
    //            os << iarr[j]<<" ";
    //        }
    //        os<<"\n";
    //    }
    return os;
}


/*
 * Print Array in a form to copy-past in PyBlastAfterglow
 */
template<class T>
static void print_xy_as_numpy(T * x_arr, double * y_arr, int nxy, int ever_i=10) {
    int i = 0;
    std::cout << "x_arr = np.array([ ";
    for (i = 0; i < nxy; i++) {
        if (i % ever_i == 0) {
            if ( (i != nxy - 1) && i != 0 )
                std::cout << ", " << x_arr[i];
            else
                std::cout << x_arr[i];
        }
    }
    std::cout << "]) \n";

    std::cout << "y_arr = np.array([ ";
    for (i = 0; i < nxy; i++) {
        if (i % ever_i == 0){
            if ( (i != nxy - 1) && i != 0 )
                std::cout << ", " << y_arr[i];
            else
                std::cout << y_arr[i];
        }
    }
    std::cout <<  "]) \n";

    std::cout << "plt.loglog(x_arr, y_arr, ls='--', label='afg')" << "\n";
    std::cout << "\n";
}
template<class T>
static void print_xy_as_numpy(T & x_arr, T & y_arr, int nxy, int ever_i=10) {
    int i = 0;
    std::cout << "x_arr = np.array([ ";
    for (i = 0; i < nxy; i++) {
        if (i % ever_i == 0) {
            if ( (i != nxy - 1) && i != 0 )
                std::cout << ", " << x_arr[i];
            else if (i > 0)
                std::cout << ", " << x_arr[i];
            else
                std::cout << x_arr[i];
        }
    }
    std::cout << "]) \n";

    std::cout << "y_arr = np.array([ ";
    for (i = 0; i < nxy; i++) {
        if (i % ever_i == 0){
            if ( (i != nxy - 1) && i != 0 )
                std::cout << ", " << y_arr[i];
            else if (i > 0)
                std::cout << ", " << y_arr[i];
            else
                std::cout << y_arr[i];
        }
    }
    std::cout <<  "]) \n";

    std::cout << "plt.loglog(x_arr, y_arr, ls='--', label='afg')" << "\n";
    std::cout << "\n";
}
template<class T>
static void print_x_as_numpy(T & x_arr, int ever_i, std::string name="x_arr", std::string label="") {
    int i = 0;
    if (!label.empty()){
        std::cout << label << "\n";
    }
    std::cout << name << " = np.array([ ";
    for (i = 0; i < x_arr.size(); i++) {
        if (i % ever_i == 0) {
            if ( (i != x_arr.size() - 1) && i != 0 )
                std::cout << ", " << x_arr[i];
            else
                std::cout << x_arr[i];
        }
    }
    std::cout << "]) \n";
}


struct Bin {
    double x_center;
    double y_center;
    double value_sum;
    int count;
};

std::vector<std::vector<Bin>> create2DHistogram(
        const std::vector<double>& x,
        const std::vector<double>& y,
        const std::vector<double>& f,
        double x_min, double x_max, int x_bins,
        double y_min, double y_max, int y_bins
) {
    // Initialize the 2D histogram bins
    std::vector<std::vector<Bin>> histogram(x_bins, std::vector<Bin>(y_bins));

    // Setup bin properties
    double x_bin_width = (x_max - x_min) / x_bins;
    double y_bin_width = (y_max - y_min) / y_bins;

    for (int i = 0; i < x_bins; ++i) {
        for (int j = 0; j < y_bins; ++j) {
            histogram[i][j].x_center = x_min + (i + 0.5) * x_bin_width;
            histogram[i][j].y_center = y_min + (j + 0.5) * y_bin_width;
            histogram[i][j].value_sum = 0.0;
            histogram[i][j].count = 0;
        }
    }

    // Fill the 2D histogram
    for (size_t i = 0; i < x.size(); ++i) {
        int x_bin = (x[i] - x_min) / x_bin_width;
        int y_bin = (y[i] - y_min) / y_bin_width;

        if (x_bin >= 0 && x_bin < x_bins && y_bin >= 0 && y_bin < y_bins) {
            histogram[x_bin][y_bin].value_sum += f[i];
            histogram[x_bin][y_bin].count++;
        }
    }

    // Now you can either leave the bins as sums, or convert them to averages by dividing with count
    // depending on your requirements.
    for (int i = 0; i < x_bins; ++i) {
        for (int j = 0; j < y_bins; ++j) {
            if (histogram[i][j].count != 0) {
                histogram[i][j].value_sum /= histogram[i][j].count;
            }
        }
    }

    return histogram;
}

std::vector<std::vector<double>> create2DHistogram2(
        const std::vector<double>& x,
        const std::vector<double>& y,
        const std::vector<double>& f,
        double x_min, double x_max, int x_bins,
        double y_min, double y_max, int y_bins
) {
    // Initialize the 2D histogram
    std::vector<std::vector<double>> histogram(x_bins, std::vector<double>(y_bins, 0.0));
    std::vector<std::vector<int>> counts(x_bins, std::vector<int>(y_bins, 0));

    double x_bin_width = (x_max - x_min) / x_bins;
    double y_bin_width = (y_max - y_min) / y_bins;

    // Fill the 2D histogram
    for (size_t i = 0; i < x.size(); ++i) {
        int x_bin = (x[i] - x_min) / x_bin_width;
        int y_bin = (y[i] - y_min) / y_bin_width;

        if ((x_bin >= 0) && (x_bin < x_bins) && (y_bin >= 0) && (y_bin < y_bins)) {
            histogram[x_bin][y_bin] += f[i];
            counts[x_bin][y_bin]++;
        }
    }

    // Convert sums to averages
    for (int i = 0; i < x_bins; ++i) {
        for (int j = 0; j < y_bins; ++j) {
            if (counts[i][j] != 0) {
                histogram[i][j] /= counts[i][j];
            }
        }
    }

    return histogram;
}

struct Output{
private:
    std::unique_ptr<logger> p_log;
public:

    static void remove_file_if_existis(const std::string &fname, std::unique_ptr<logger> & p_log){

//    LOG_INIT_CERR();
//    Log.set_log_level(LOG_SILENT);

        try {
            if (std::experimental::filesystem::remove(fname))
//            (*p_log)(LOG_INFO) << "file " << fname << " deleted.\n";
                (*p_log)(LOG_INFO,AT) << "file " << fname << " deleted.\n";
            else
//            (*p_log)(LOG_INFO) << "file " << fname << " not found.\n";
                (*p_log)(LOG_INFO,AT) << "file " << fname << " not found.\n";
        }
        catch(const std::experimental::filesystem::filesystem_error& err) {
//        (*p_log)(LOG_ERR) << "filesystem error: " << err.what() << '\n';
            (*p_log)(LOG_ERR,AT) << " filesystem error: " << err.what() << '\n';
//        std::cerr <<AT<<"\n";
        }
    }

    explicit Output(int loglevel){
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "Output");
    }

    template<typename vecarr>
    void Data(vecarr & data,
              std::vector<std::string> &names,
              const std::string& fpath,
              const bool show_output = true){

//    LOG_INIT_CERR();
//    Log.set_log_level(LOG_SILENT);

        // remove the old file (so to avoid data piling up)
        remove_file_if_existis(fpath, p_log);

        std::fstream of(fpath, std::ios::out | std::ios::app);

        if (of.is_open())
        {
            stream_arr(of, data, names);
            if (show_output)
                stream_arr(std::cout, data, names);
            of.close();
        } else {
            std::cout << "Error! Could not open the file " << fpath << "\n";
        }

    }


    /**
 * Save vector of 2D vectors in a single h5
 * WARNING! using h5Cpp.h library!
 * @param data
 * @param set_names
 * @param fpath
 */
    void VectorOfTablesH5(std::vector<VecVector> & data,
                          std::vector<std::string> set_names,
                          std::string fpath){

//    std::cout << "S"

        // remove the old file (so to avoid data piling up)
        remove_file_if_existis(fpath, p_log);

        H5::H5File file(fpath, H5F_ACC_TRUNC); // "/home/m/Desktop/tryout/file.h5"

        size_t n_entries = data.size();

        for (size_t id = 0; id < n_entries; ++id){

            size_t nrow = data[id].size();
            size_t ncol = data[id][0].size();

            /// convert VecVector to array[][]
            double varray[nrow][ncol];
            for( size_t i = 0; i<nrow; ++i) {
                for( size_t j = 0; j<ncol; ++j) {
                    varray[i][j] = data[id][i][j];
                }
            }

            /// create dimensions
            hsize_t dimsf[2];
            dimsf[0] = nrow;
            dimsf[1] = ncol;
            H5::DataSpace dataspace(2, dimsf);
            H5::DataType datatype(H5::PredType::NATIVE_DOUBLE);
            H5::DataSet dataset = file.createDataSet(set_names[id], datatype, dataspace);
            dataset.write(varray, H5::PredType::NATIVE_DOUBLE);
            dataset.close();
            dataspace.close();
        }
        file.close();
    }


    static void addVector(Vector & data, std::string name, H5::H5File & file){
        if(data.empty()){ std::cerr << AT << " no data\n"; exit(1); }
        size_t nrow = data.size();

        /// convert VecVector to array[][]
        double varray[nrow];
        for( size_t i = 0; i<nrow; ++i) { varray[i] = data[i]; }

        // preparation of a dataset and a file.
        hsize_t dimsf[1] = {nrow };
        H5::DataSpace dataspace(1, dimsf);
        H5::DataType datatype(H5::PredType::NATIVE_DOUBLE);
        H5::DataSet dataset = file.createDataSet(name, datatype, dataspace);
        dataset.write( varray, H5::PredType::NATIVE_DOUBLE);
        dataset.close();
        dataspace.close();
    }

    static void add2Dtable(VecVector & data, std::string name, H5::H5File & file){
        if(data.empty()){ std::cerr << AT << " no data\n"; exit(1); }
        if(data[0].empty()){ std::cerr << AT << " no data[0]\n"; exit(1); }

        size_t nrow = data.size();
        size_t ncol = data[0].size();

        /// convert VecVector to array[][]
        double varray[nrow][ncol];
        for( size_t i = 0; i<nrow; ++i) {
            for( size_t j = 0; j<ncol; ++j) {
                varray[i][j] = data[i][j];
            }
        }

        /// create dimensions
        hsize_t dimsf[2];
        dimsf[0] = nrow;
        dimsf[1] = ncol;
        H5::DataSpace dataspace(2, dimsf);
        H5::DataType datatype(H5::PredType::NATIVE_DOUBLE);
        H5::DataSet dataset = file.createDataSet(name, datatype, dataspace);
        dataset.write(varray, H5::PredType::NATIVE_DOUBLE);
        dataset.close();
        dataspace.close();
    }

    static void addStrDbMap(StrDbMap & attrs, H5::H5File & file){
        for (auto & ele : attrs){
            auto key = ele.first;
            auto val = ele.second;
            H5::IntType int_type(H5::PredType::NATIVE_DOUBLE);
            H5::DataSpace att_space(H5S_SCALAR);
            H5::Attribute att = file.createAttribute(key, int_type, att_space );
            att.write( int_type, &val );
        }
    }

    static void addGroupWith1Ddata(VecVector & data, std::string group_name, std::vector<std::string> array_names, H5::H5File & file){
        if(data.empty()){ std::cerr << AT << " no data\n"; exit(1); }
        if(data[0].empty()){ std::cerr << AT << " no data[0]\n"; exit(1); }
        if(array_names.empty()){ std::cerr << AT << " no array_names\n"; exit(1); }
        if(data.size()!=array_names.size()){
            std::cerr << data.size() << ": " << array_names.size()<<"\n";
            std::cerr << AT << " data.size()!=array_names.size()\n";
            exit(1);
        }
        H5::Group grp(file.createGroup(group_name));
        for (size_t iv = 0; iv < data.size(); iv++){
            double varray[data[iv].size()];
            for (size_t ii = 0; ii < data[iv].size(); ii++)
                varray[ii] = data[iv][ii];
            hsize_t dimsf[1];
            dimsf[0] = data[iv].size();
            H5::DataSpace dataspace(1, dimsf);
            H5::DataType datatype(H5::PredType::NATIVE_DOUBLE);
            H5::DataSet dataset = grp.createDataSet(array_names[iv], datatype, dataspace);
            dataset.write(varray, H5::PredType::NATIVE_DOUBLE);
            dataset.close();
            dataspace.close();
        }
        grp.close();
    }

    static void addGroupWith2Ddata(std::vector<VecVector> & data, std::string group_name, std::vector<std::string> table_names, H5::H5File & file){
        if(data.empty()){ std::cerr << AT << " no data\n"; exit(1); }
        if(data[0].empty()){ std::cerr << AT << " no data[0]\n"; exit(1); }
        if(table_names.empty()){ std::cerr << AT << " no array_names\n"; exit(1); }
        if(data.size()!=table_names.size()){ std::cerr << AT << " data.size()!=table_names.size()\n"; exit(1); }

        H5::Group grp(file.createGroup(group_name));
        for (size_t iv = 0; iv < data.size(); iv++){
            size_t nrow = data[iv].size();
            size_t ncol = data[iv][0].size(); /// ncoll is now different for each row

            auto * varray = new double[nrow*ncol];

            size_t kk = 0;
            for (size_t ii = 0; ii < nrow; ii++)
                for(size_t jj = 0; jj < ncol; jj++) {
                    varray[kk] = data[iv][ii][jj];
                    kk++;
                }

            hsize_t dimsf[2];
            dimsf[0] = nrow;
            dimsf[1] = ncol;
            H5::DataSpace dataspace(2, dimsf);
            H5::DataType datatype(H5::PredType::NATIVE_DOUBLE);
            H5::DataSet dataset = grp.createDataSet(table_names[iv], datatype, dataspace);
            dataset.write(varray, H5::PredType::NATIVE_DOUBLE);
            dataset.close();
            dataspace.close();

            delete[] varray;
        }
        grp.close();
    }







    void VectorOfVectorsH5(VecVector & data,
                           std::vector<std::string> set_names,
                           std::string fpath,
                           std::unordered_map<std::string, double> attrs={}){

//    std::cout << "S"
        if (data.empty()){
            (*p_log)(LOG_ERR,AT) <<" Empty data cannot be saved.\n"
                      << " Exiting...\n";
//            std::cerr << AT << "\n";
            exit(1);
        }
        if (data[0].empty()){
            (*p_log)(LOG_ERR,AT) <<" Empty data[0] cannot be saved\n"
                      << " Exiting...\n";
//            std::cerr << AT << "\n";
            exit(1);
        }
        if (data.size() != set_names.size()){
            (*p_log)(LOG_ERR,AT) <<"Datasize and name.size are not the same; Data has "<<data.size()<<" names have "<<set_names.size()<<"\n"
                      << " Exiting...\n";
//            std::cerr << AT << "\n";
            exit(1);
        }

        // remove the old file (so to avoid data piling up)
        remove_file_if_existis(fpath, p_log);

        H5::H5File file(fpath, H5F_ACC_TRUNC); // "/home/m/Desktop/tryout/file.h5"

        size_t n_entries = data.size();

        for (size_t id = 0; id < n_entries; ++id){

            size_t nrow = data[id].size();

            /// convert VecVector to array[][]
            double varray[nrow];
            for( size_t i = 0; i<nrow; ++i) { varray[i] = data[id][i]; }

            // preparation of a dataset and a file.
            hsize_t dimsf[1] = {nrow };
            H5::DataSpace dataspace(1, dimsf);
            H5::DataType datatype(H5::PredType::NATIVE_DOUBLE);
            H5::DataSet dataset = file.createDataSet(set_names[id], datatype, dataspace);
            dataset.write( varray, H5::PredType::NATIVE_DOUBLE);
            dataset.close();
            dataspace.close();
        }

        if (!attrs.empty()){
            for (auto & ele : attrs){
                auto key = ele.first;
                auto val = ele.second;
                H5::IntType int_type(H5::PredType::NATIVE_DOUBLE);
                H5::DataSpace att_space(H5S_SCALAR);
                H5::Attribute att = file.createAttribute(key, int_type, att_space );
                att.write( int_type, &val );
            }
        }

        file.close();
    }

    void VecVectorOfVectorsAsGroupsH5(std::vector<VecVector>  & data,
                                      std::vector<std::string> group_names,
                                      std::vector<std::string> arr_names,
                                      std::string fpath,
                                      std::unordered_map<std::string, double> attrs = {},
                                      std::vector<std::unordered_map<std::string, double>> group_attrs = {}
    ){
        // remove the old file (so to avoid data piling up)
        remove_file_if_existis(fpath, p_log);

        H5::H5File file(fpath, H5F_ACC_TRUNC); // "/home/m/Desktop/tryout/file.h5"

        size_t n_entries = data.size();

        for (size_t id = 0; id < n_entries; ++id){
            H5::Group grp(file.createGroup(group_names[id]));
            for (size_t iv = 0; iv < data[id].size(); iv++){
                double varray[data[id][iv].size()];
                for (size_t ii = 0; ii < data[id][iv].size(); ii++)
                    varray[ii] = data[id][iv][ii];
                hsize_t dimsf[1];
                dimsf[0] = data[id][iv].size();
                H5::DataSpace dataspace(1, dimsf);
                H5::DataType datatype(H5::PredType::NATIVE_DOUBLE);
                H5::DataSet dataset = grp.createDataSet(arr_names[iv], datatype, dataspace);
                dataset.write(varray, H5::PredType::NATIVE_DOUBLE);
                dataset.close();
                dataspace.close();
            }
            // write attributes for the data
            if (!group_attrs.empty()){
                if (group_attrs.size() != group_names.size()){
                    (*p_log)(LOG_ERR,AT)  << " size mismatch group_attrs="<<group_attrs.size()<<" group_names="<<group_names.size()<<"\n"
                               << " Exiting...\n";
//                    std::cerr << AT << "\n";
                    exit(1);
                }
                auto & group_attr = group_attrs[id];
                for (auto & ele : group_attr){
                    auto key = ele.first;
                    auto val = ele.second;
                    H5::IntType int_type(H5::PredType::NATIVE_DOUBLE);
                    H5::DataSpace att_space(H5S_SCALAR);
                    H5::Attribute att = grp.createAttribute(key, int_type, att_space );
                    att.write( int_type, &val );
                }
            }
            grp.close();
        }

        // saving attrs

        if (!attrs.empty()){
            for (auto & ele : attrs){
                auto key = ele.first;
                auto val = ele.second;
                H5::IntType int_type(H5::PredType::NATIVE_DOUBLE);
                H5::DataSpace att_space(H5S_SCALAR);
                H5::Attribute att = file.createAttribute(key, int_type, att_space );
                att.write( int_type, &val );
            }
        }

        file.close();
    }

    void VecVectorOfVectorsAsGroupsAndVectorOfVectorsH5(
            std::string fpath,
            std::vector<VecVector> & group_data,
            std::vector<std::string> group_names,
            std::vector<std::string> vector_names_in_each_group,
            std::vector<std::vector<double>> & data,
            std::vector<std::string> vector_names,
            std::unordered_map<std::string, double> attrs = {}
    ){
        // remove the old file (so to avoid group_data piling up)
        remove_file_if_existis(fpath, p_log);

        H5::H5File file(fpath, H5F_ACC_TRUNC); // "/home/m/Desktop/tryout/file.h5"

        /// save group data

        size_t n_entries = group_data.size();

        for (size_t id = 0; id < n_entries; ++id){
            H5::Group grp(file.createGroup(group_names[id]));
            for (size_t iv = 0; iv < group_data[id].size(); iv++){
                double varray[group_data[id][iv].size()];
                for (size_t ii = 0; ii < group_data[id][iv].size(); ii++)
                    varray[ii] = group_data[id][iv][ii];
                hsize_t dimsf[1];
                dimsf[0] = group_data[id][iv].size();
                H5::DataSpace dataspace(1, dimsf);
                H5::DataType datatype(H5::PredType::NATIVE_DOUBLE);
                H5::DataSet dataset = grp.createDataSet(vector_names_in_each_group[iv], datatype, dataspace);
                dataset.write(varray, H5::PredType::NATIVE_DOUBLE);
                dataset.close();
                dataspace.close();
            }
            grp.close();
        }

        /// save separate data

        if (data.empty()){
            (*p_log)(LOG_ERR,AT) <<" Empty data cannot be saved. "
                      << " Exiting...\n";
//            std::cerr << AT << "\n";
            exit(1);
        }
        if (data[0].empty()){
            (*p_log)(LOG_ERR,AT) <<" Empty data[0] cannot be saved\n"
                      << " Exiting...\n";
//            std::cerr << AT << "\n";
            exit(1);
        }
        if (data.size() != vector_names.size()){
            (*p_log)(LOG_ERR,AT) <<"Datasize and name.size are not the same; Data has "
                      <<data.size()<<" names have "<<vector_names.size()
                      << " Exiting...\n";
//            std::cerr << AT << "\n";
            exit(1);
        }

        size_t n_entries_ = data.size();

        for (size_t id = 0; id < n_entries_; ++id){

            size_t nrow = data[id].size();
//        size_t ncol = data[id][0].size();

            /// convert VecVector to array[][]
            double varray[nrow];
            for( size_t i = 0; i<nrow; ++i) { varray[i] = data[id][i]; }

            // preparation of a dataset and a file.
//        hsize_t dim[1];
//        dim[0] = data.size();                   // using vector::size()
//        int rank = sizeof(dim) / sizeof(hsize_t);
//        H5::DataSpace dataspace(rank, dim);
//
            hsize_t dimsf[1] = {nrow };
//        dimsf[0] = nrow;
//        dimsf[1] = ncol;
            H5::DataSpace dataspace(1, dimsf);
            H5::DataType datatype(H5::PredType::NATIVE_DOUBLE);
            H5::DataSet dataset = file.createDataSet(vector_names[id], datatype, dataspace);
            dataset.write( varray, H5::PredType::NATIVE_DOUBLE);
            dataset.close();
            dataspace.close();
        }

        // saving attrs

        if (!attrs.empty()){
            for (auto & ele : attrs){
                auto key = ele.first;
                auto val = ele.second;
                H5::IntType int_type(H5::PredType::NATIVE_DOUBLE);
                H5::DataSpace att_space(H5S_SCALAR);
                H5::Attribute att = file.createAttribute(key, int_type, att_space );
                att.write( int_type, &val );
            }
        }

        file.close();
    }

    void VectorOfTablesAsGroupsAndVectorOfVectorsH5(
            std::string fpath,
            std::vector<std::vector<VecVector>> & group_data,
            std::vector<std::string> group_names_tables,
            std::vector<std::string> table_names_in_each_group,
            std::vector<std::vector<double>> & data,
            std::vector<std::string> vector_names = {},
            std::unordered_map<std::string, double> attrs = {}
    ){
        // remove the old file (so to avoid group_data piling up)
        remove_file_if_existis(fpath, p_log);

        H5::H5File file(fpath, H5F_ACC_TRUNC); // "/home/m/Desktop/tryout/file.h5"

        /// save group data
        if (group_names_tables.size()!=group_data.size()){
            (*p_log)(LOG_ERR,AT)  << " Group names="<<group_names_tables.size()
                       <<" while data vec contains=" << group_data.size() << " groups of tables"
                       <<" Exiting...";
//            std::cerr << AT << "\n";
            exit(1);
        }
        if (table_names_in_each_group.size()!=group_data[0].size()){
            (*p_log)(LOG_ERR,AT)  << " Table names="<<table_names_in_each_group.size()
                       <<" while each group contains=" << group_data[0].size() << " tables"
                       << " Exiting...\n";
//            std::cerr << AT << "\n";
            exit(1);
        }

        size_t n_entries = group_data.size();

//    std::cout << " adding tables groups to h5\n";

        // allocate buffer storage

//    double ***varr_ = new double**[n_entries];
//    for (size_t iv = 0; n_entries; iv++) {
//        varr_[iv] = new double*[group_data[iv].size()];
//        for (j = 0; j < group_data[id][iv][0].size(); ++j)
//            arr[i][j] = new int[Z];
//    }

        /// /// output container
        //        std::vector< // times & freqs
        //                std::vector< // v_ns
        //                        std::vector< // shells
        //                                std::vector<double>>>> // data
        //        out {};

        for (size_t id = 0; id < n_entries; ++id){
            H5::Group grp(file.createGroup(group_names_tables[id]));
            for (size_t iv = 0; iv < group_data[id].size(); iv++){
                size_t nrow = group_data[id][iv].size();
                size_t ncol = group_data[id][iv][0].size(); /// ncoll is now different for each row

//            auto** varray = new double*[nrow];
//            for(int i = 0; i < nrow; ++i)
//                varray[i] = new double[ncol];

//            double (*varray)[ncols] = malloc(nrow * sizeof(*varray));
//            double *varray = new double[nrow*ncol];


//            double varray[nrow*ncol];
                auto * varray = new double[nrow*ncol];

                size_t kk = 0;
                for (size_t ii = 0; ii < nrow; ii++)
                    for(size_t jj = 0; jj < ncol; jj++) {
//                    varray[ii] = group_data[id][iv][ii][jj];
                        varray[kk] = group_data[id][iv][ii][jj];
                        kk++;
                    }
//            size_t nrow = data[id].size();
//            size_t ncol = data[id][0].size();
//
//            /// convert VecVector to array[][]
//            double varray[nrow][ncol];
//            for( size_t i = 0; i<nrow; ++i) {
//                for( size_t j = 0; j<ncol; ++j) {
//                    varray[i][j] = data[id][i][j];
//                }
//            }
//            double *outdata = (double *)calloc(nrow*ncol, sizeof(double));
//            for (size_t ii = 0; ii < nrow; ii++)
//                for(size_t jj = 0; jj < ncol; jj++) {
//                    *(data + i * ncol + j) = (ii * ncol + jj);
//                    outdata[ii][jj] = group_data[id][iv][ii][jj];
//                }




                /// create dimensions
//            std::cout << " saving.." << table_names_in_each_group[iv] << "\n";
//            std::cout << varray[0][0] << "\n";
//            std::cout << varray[0][10] << "\n";
//            std::cout << varray[0][100] << "\n";
//            std::cout << varray[0][1000] << "\n";
                hsize_t dimsf[2];
                dimsf[0] = nrow;
                dimsf[1] = ncol;
                H5::DataSpace dataspace(2, dimsf);
                H5::DataType datatype(H5::PredType::NATIVE_DOUBLE);
                H5::DataSet dataset = grp.createDataSet(table_names_in_each_group[iv], datatype, dataspace);
                dataset.write(varray, H5::PredType::NATIVE_DOUBLE);
                dataset.close();
                dataspace.close();

                delete[] varray;


//            hsize_t dimsf[1];
//            dimsf[0] = group_data[id][iv].size();
//            H5::DataSpace dataspace(1, dimsf);
//            H5::DataType datatype(H5::PredType::NATIVE_DOUBLE);
//            H5::DataSet dataset = grp.createDataSet(vector_names_in_each_group[iv], datatype, dataspace);
//            dataset.write(varray, H5::PredType::NATIVE_DOUBLE);
//            dataset.close();
//            dataspace.close();

            }
            grp.close();


        }

        /// save separate data
        if (data.empty()) {
            (*p_log)(LOG_ERR,AT)  << " Empty data cannot be saved.\n";
//            std::cerr << AT << "\n";
            exit(1);
        }
        if (data[0].empty()) {
            (*p_log)(LOG_ERR,AT)  << " Empty data[0] cannot be saved\n";
//            std::cerr << AT << "\n";
            exit(1);
        }
        if (data.size() != vector_names.size()) {
            (*p_log)(LOG_ERR,AT)  << "Datasize and name.size are not the same; Data has " << data.size() << " names have "
                       << vector_names.size() << "\n";
//            std::cerr << AT << "\n";
            exit(1);
        }

        size_t n_entries_ = data.size();

//    std::cout << " adding other arrays to h5\n";
        for (size_t id = 0; id < n_entries_; ++id) {

            size_t nrow = data[id].size();
//        size_t ncol = data[id][0].size();

            /// convert VecVector to array[][]
            double varray[nrow];
            for (size_t i = 0; i < nrow; ++i) { varray[i] = data[id][i]; }

            // preparation of a dataset and a file.
//        hsize_t dim[1];
//        dim[0] = data.size();                   // using vector::size()
//        int rank = sizeof(dim) / sizeof(hsize_t);
//        H5::DataSpace dataspace(rank, dim);
//
            hsize_t dimsf[1] = {nrow};
//        dimsf[0] = nrow;
//        dimsf[1] = ncol;
            H5::DataSpace dataspace(1, dimsf);
            H5::DataType datatype(H5::PredType::NATIVE_DOUBLE);
            auto _name = vector_names[id];
            H5::DataSet dataset = file.createDataSet(_name, datatype, dataspace);
            dataset.write(varray, H5::PredType::NATIVE_DOUBLE);
            dataset.close();
            dataspace.close();
        }

//    std::cout << " adding attrs to h5\n";

        // saving attrs
        if (!attrs.empty()){
            for (auto & ele : attrs){
                auto key = ele.first;
                auto val = ele.second;
                H5::IntType int_type(H5::PredType::NATIVE_DOUBLE);
                H5::DataSpace att_space(H5S_SCALAR);
                H5::Attribute att = file.createAttribute(key, int_type, att_space );
                att.write( int_type, &val );
            }
        }
//    std::cout << " closing file  h5\n";
        file.close();
//    std::cout << " done\n";
    }

};

/// https://www.techiedelight.com/print-keys-values-map-cpp/
class myclass
{
    template<typename K, typename V>
    void operator()(const std::pair<K, V> &p) {
        std::cout << "{" << p.first << ": " << p.second << "}\n";
    }
} ob;

template<typename K, typename V>
void print(const std::pair<K, V> &p) {
    std::cout << "{" << p.first << ": " << p.second << "}\n";
}

template<typename K, typename V>
void print_map(std::unordered_map<K, V> const &m)
{
    // specify a lambda expression
    std::for_each(m.begin(),
                  m.end(),
                  [](const std::pair<int, char> &p) {
                      std::cout << "{" << p.first << ": " << p.second << "}\n";
                  });

    // or pass an object of a class overloading the ()operator
    // std::for_each(m.begin(), m.end(), ob);

    // or specify a function
    // std::for_each(m.begin(), m.end(), print<K, V>);
}

void cast_times_freqs(Vector& lc_times, Vector& lc_freqs,
                      Vector& _times, Vector& _freqs,
                      bool is_one_to_one_already, std::unique_ptr<logger> & p_log){
    if (lc_times.empty() || lc_freqs.empty()){
        (*p_log)(LOG_ERR,AT)<<" isEmpty time or freq arr.\n";
        exit(1);
    }
    if (is_one_to_one_already){
        if (lc_times.size()!=lc_freqs.size()){
            (*p_log)(LOG_ERR,AT)<<" size mismatch between arrays time and freq (for one-to-one freq-to-time)\n";
            exit(1);
        }

        _times = lc_times;
        _freqs = lc_freqs;

    }
    else {
        _times.resize(lc_freqs.size() * lc_times.size(), 0.0);
        _freqs.resize(lc_freqs.size() * lc_times.size(), 0.0);
        size_t ii = 0;
        for (double time: lc_times) {
            for (double freq: lc_freqs) {
                _times[ii] = time;
                _freqs[ii] = freq;
                ii++;
            }
        }
    }
}

void readParFile2(std::unordered_map<std::string, double> & pars,
                  std::unordered_map<std::string, std::string> & opts,
                  std::unique_ptr<logger> & p_log,
                  std::string parfile_path,
                  std::string from_line, std::string until_line
){

    /// settings for reading the parfile
    std::string key_after_which_to_look_for_parameters = "* Parameters";
    std::string key_after_which_to_look_for_settings = "* Settings";
    char char_that_separaters_name_and_value = '=';
    char char_that_separaters_value_and_comment = '#';
    std::vector<std::string> leave_spaces_for = {
            "lc_freqs", "lc_times", "skymap_freqs", "skymap_times", "spec_times", "spec_freqs", "spec_gams"
    };

//    std::unique_ptr<logger> p_log;
//    p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "readParFile2");
    if (!std::experimental::filesystem::exists(parfile_path)) {
        (*p_log)(LOG_ERR, AT) << " Parfile not found. " + parfile_path << "\n";
        exit(1);

    }
    std::ifstream fin(parfile_path);
    std::string line;

    bool is_in_the_reqired_block = false;
    bool reading_pars = false;
    bool reading_opts = false;
    while (std::getline(fin, line)) {
        /// check if reading the required block of parfile
        if (line == from_line)
            is_in_the_reqired_block = true;
        if (line == until_line)
            is_in_the_reqired_block = false;
        if (!is_in_the_reqired_block)
            continue;

        /// read parameters (double) and settings (str) separately
        if (line == key_after_which_to_look_for_parameters) {
            reading_pars = true; reading_opts = false;
        }
        if (line == key_after_which_to_look_for_settings) {
            reading_pars = false; reading_opts = true;
        }
        if (line[0] == char_that_separaters_value_and_comment)
            continue;

        /// read pars (str, double)
        if (reading_pars and (line.length() > 1) and (line != key_after_which_to_look_for_parameters)) {
            unsigned long pos = line.find_first_of(char_that_separaters_name_and_value);
            std::string val = line.substr(pos + 1), par = line.substr(0, pos);
            par.erase(std::remove_if(par.begin(), par.end(), ::isspace), par.end());
//            val.erase(std::remove_if(val.begin(), val.end(), ::isspace), val.end());
            if (val.find(char_that_separaters_value_and_comment) != std::string::npos){
                unsigned long _pos = val.find_first_of(char_that_separaters_value_and_comment);
                std::string _comment = val.substr(_pos + 1), _val = val.substr(0, _pos);
                val = _val;
            }
            val.erase(std::remove_if(val.begin(), val.end(), ::isspace), val.end());
            double value = std::stod(val);
            pars.insert(std::pair<std::string, double>(par, value));
        }
        /// read opts (str, str)
        if (reading_opts and (line.length() > 1) and (line != key_after_which_to_look_for_settings)) {
            unsigned long pos = line.find_first_of(char_that_separaters_name_and_value);
            std::string val = line.substr(pos + 1),
                    par = line.substr(0, pos);
            par.erase(std::remove_if(par.begin(), par.end(), ::isspace), par.end());
//            val.erase(std::remove_if(val.begin(), val.end(), ::isspace), val.end());
            if (val.find(char_that_separaters_value_and_comment) != std::string::npos){
                unsigned long _pos = val.find_first_of(char_that_separaters_value_and_comment);
                std::string _comment = val.substr(_pos + 1), _val = val.substr(0, _pos);
                val = _val;
            }
            if (std::find(leave_spaces_for.begin(), leave_spaces_for.end(), par) == leave_spaces_for.end())
                val.erase(std::remove_if(val.begin(), val.end(), ::isspace), val.end());
//            double value = std::stod(val);
//            if (val[0]=='a'&&val[1]=='r'&&val[2]=='r'&&val[3]=='a'&&val[4]=='y'){
//                std::cout << par << " = " << val<<"\n";
//            }
            opts.insert(std::pair<std::string, std::string>(par, val));
        }
    }

}

Vector makeVecFromString(const std::string line, std::unique_ptr<logger> & p_log){
//    std::string line = _line;
    char space_char = ' ';
//    char end_char = '#';
    std::vector<std::string> words{};
    std::stringstream sstream(line);
    std::string word;
    while (std::getline(sstream, word, space_char)){
//        word.erase(std::remove_if(word.begin(), word.end(), ispunct), word.end()); // TODO removes '.' from floats
        words.push_back(word);
    }
    if (words.size() == 1){
        (*p_log)(LOG_ERR, AT) << "incorrect word size in array line. Expected 2 or 3 found only " << words[0] << '\n';
        exit(1);
    }
    /// remove first emtpy element left after '= ' this
    if ((words[0].empty())or(words[0] == " ")){
        words.erase(words.begin());
    }

    /// construct the vector
    std::vector<std::string> words2{};
    words2.reserve(words.size());
    for (auto _word : words) words2.emplace_back(_word);

    Vector res;
    if (words[0] == "array"){
        if (words[1] == "logspace"){
            const double val1 = std::log10( std::stod(words[2]) );
            const double val2 = std::log10( std::stod(words[3]) );
            const int nvals = (int)std::stod(words[4]);
            res = TOOLS::MakeLogspaceVec(val1,val2,nvals,10);
        }
        else{
            words.erase(words.begin());
            for (const auto &str : words){
                res.push_back(std::stod(str));
            }
        }
    }
    else{
        (*p_log)(LOG_ERR, AT) << "incorrect first word in array line. Expected 'array' found " << words[0] << '\n';
        exit(1);
    }
//    std::cout << res << "\n";
    return std::move( res );
}

std::vector<int> findDuplicateCoordinates(const std::vector<double>& x_coords, const std::vector<double>& y_coords) {
    if (x_coords.size() != y_coords.size()) {
        throw std::runtime_error("Mismatched vector sizes.");
    }

    // Pair up the coordinates with their original indices.
    std::vector<std::pair<std::pair<double, double>, int>> coordsWithIndices;
    for (size_t i = 0; i < x_coords.size(); ++i) {
        coordsWithIndices.push_back({{x_coords[i], y_coords[i]}, static_cast<int>(i)});
    }

    // Sort by the coordinates.
    std::sort(coordsWithIndices.begin(), coordsWithIndices.end());

    // Find duplicates
    std::vector<int> duplicateIndices;
    for (size_t i = 1; i < coordsWithIndices.size(); ++i) {
        if (coordsWithIndices[i].first == coordsWithIndices[i - 1].first) {
            duplicateIndices.push_back(coordsWithIndices[i].second);
            if (i == coordsWithIndices.size() - 1 || coordsWithIndices[i].first != coordsWithIndices[i + 1].first) {
                duplicateIndices.push_back(coordsWithIndices[i - 1].second);
            }
        }
    }

    // Sort indices to return them in the original order.
    std::sort(duplicateIndices.begin(), duplicateIndices.end());

    return duplicateIndices;
}

/// https://www.geeksforgeeks.org/LinearRegression-analysis-and-the-best-fitting-line-using-c/
class LinearRegression {
    // Dynamic array which is going
    // to contain all (i-th x)
    Vector & x;

    // Dynamic array which is going
    // to contain all (i-th y)
    Vector & y;

    // final non-zero data index
    // in case x,y have max indx.
    size_t ir=0;

    // Store the coefficient/slope in
    // the best fitting line
    double coeff{};

    // Store the constant term in
    // the best fitting line
    double constTerm{};

    // Contains sum of product of
    // all (i-th x) and (i-th y)
    double sum_xy{};

    // Contains sum of all (i-th x)
    double sum_x{};

    // Contains sum of all (i-th y)
    double sum_y{};

    // Contains sum of square of
    // all (i-th x)
    double sum_x_square{};

    // Contains sum of square of
    // all (i-th y)
    double sum_y_square{};

public:
    // Constructor to provide the default
    // values to all the terms in the
    // object of class LinearRegression
    LinearRegression(Vector & x, Vector & y) : x(x), y(y){
        coeff = 0;
        constTerm = 0;
        sum_y = 0;
        sum_y_square = 0;
        sum_x_square = 0;
        sum_x = 0;
        sum_xy = 0;
    }

    // Function that calculate the coefficient/
    // slope of the best fitting line
    void calculateCoefficient()
    {
//        auto N = (double)x.size();
        double numerator = ((double)ir * sum_xy - sum_x * sum_y);
        double denominator = ((double)ir * sum_x_square - sum_x * sum_x);
        coeff = numerator / denominator;
    }

    // Member function that will calculate
    // the constant term of the best
    // fitting line
    void calculateConstantTerm()
    {
//        auto N = (double)x.size();
        double numerator = (sum_y * sum_x_square - sum_x * sum_xy);
        double denominator = ((double)ir * sum_x_square - sum_x * sum_x);
        constTerm = numerator / denominator;
    }

    // Function that return the number
    // of entries (xi, yi) in the data set
//    size_t sizeOfData() {
//        return x.size();
//    }

    // Function that return the coefficient/
    // slope of the best fitting line
    double coefficient() {
        if (coeff == 0)
            calculateCoefficient();
        return coeff;
    }
    double getCoeff(){return coeff;}

    // Function that return the constant
    // term of the best fitting line
    double constant() {
        if (constTerm == 0)
            calculateConstantTerm();
        return constTerm;
    }
    double getConst(){return constTerm;}

    bool isTrained(){
        if ((coeff == 0) &&( constTerm == 0)) return false;
        return true;
    }

    // Function that print the best
    // fitting line
    void PrintBestFittingLine() {
        if (coeff == 0 && constTerm == 0) {
            calculateCoefficient();
            calculateConstantTerm();
        }
        std::cout << "The best fitting line is y = "
             << coeff << "x + " << constTerm << std::endl;
    }

    // Function to take input from the dataset
    void commputeInput() {
        for (size_t i = 0; i < x.size(); i++) {
            if ((x[i] != 0.) and (y[i] != 0.)) {
                sum_xy += x[i] * y[i];
                sum_x += x[i];
                sum_y += y[i];
                sum_x_square += x[i] * x[i];
                sum_y_square += y[i] * y[i];
                ir++;
            }
        }
    }

    // Function to show the data set
    void showData()
    {
        for (size_t i = 0; i < 62; i++) {
            printf("_");
        }
        printf("\n\n");
        printf("|%15s%5s %15s%5s%20s\n",
               "X", "", "Y", "", "|");

        for (size_t i = 0; i < x.size(); i++) {
            printf("|%20f %20f%20s\n",
                   x[i], y[i], "|");
        }

        for (size_t i = 0; i < 62; i++) {
            printf("_");
        }
        printf("\n");
    }

    // Function to predict the value
    // corresponding to some input
    double predict(double x) const {
        return coeff * x + constTerm;
    }

    // Function that returns overall
    // sum of square of errors
    double errorSquare() {
        double ans = 0;
        for (size_t i = 0; i < ir; i++) {
            ans += ((predict(x[i]) - y[i])
                    * (predict(x[i]) - y[i]));
        }
        return ans;
    }

    // Functions that return the error
    // i.e the difference between the
    // actual value and value predicted
    // by our model
    double errorIn(double num) {
        for (size_t i = 0; i < ir; i++) {
            if (num == x[i]) {
                return (y[i] - predict(x[i]));
            }
        }
        return 0;
    }
};

#endif //SRC_UTILS_H
