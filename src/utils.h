//
// Created by vsevolod on 21/12/22.
//

#ifndef SRC_UTILS_H
#define SRC_UTILS_H

#include "pch.h"
#include "logger.h"

typedef std::valarray<double> Array;
typedef std::vector<Array> VecArray;
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
    static std::valarray<double> MakeLogspace(const double &start,
                                              const double &stop,
                                              const int &num = 50,
                                              const double &base = 10) {
        auto retval = MakeLogspaceVec(start, stop, num, base);
        return std::move(std::valarray<double> (retval.data(), retval.size()));

    }

    std::vector<double> linspace(double first, double last, int len) {
        std::vector<double> result(len);
        double step = (last-first) / (len - 1);
        for (int i=0; i<len; i++) { result[i] = first + i*step; }
        return result;
    }

}

bool getBoolOpt(std::string opt, std::unordered_map<std::string,std::string> & opts,
                const char *loc, std::unique_ptr<logger> & p_log, bool def = false, bool req = false){
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
std::string getStrOpt(std::string opt, std::unordered_map<std::string,std::string> & opts,
                      const char *loc, std::unique_ptr<logger> & p_log, std::string def = "", bool req = false){
//    opt = "use_dens_prof_behind_jet_for_ejecta";
    std::string b_val;
    std::string val;
    if ( opts.find(opt) == opts.end() )
        if (req){
            (*p_log)(LOG_ERR, loc) << "Required parameter " << opt << " not found. exiting...\n";
//            std::cerr << AT << " \n Required parameter " << par << " not found. exiting...\n";
            exit(1);
        }
        else{
            b_val = def;
        }
    else{
        val = (std::string)opts.at(opt);
    }
    return val;
}
double getDoublePar(std::string par, std::unordered_map<std::string,double> & pars,
                    const char *loc, std::unique_ptr<logger> & p_log, double def = 0., bool req = false){
    double val;
    if ( pars.find(par) == pars.end() ){
        if (req){
            (*p_log)(LOG_ERR, loc) << " Required parameter " << par << " not found. exiting...\n";
            std::cerr << " Required parameter " << par << " not found. exiting...\n";
            exit(1);
        }
        else{
            val = def;
        }
    }
    else{
        val = (double)pars.at(par);
    }
    return val;
}

Vector arrToVec(Array & array){
    std::vector<double> vec;
    vec.assign(std::begin(array), std::end(array));
    return std::move( vec );
}

/// https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {

    // initialize original index locations
    std::vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values
    stable_sort(idx.begin(), idx.end(),
                [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

    return std::move( idx );
}
template <typename T>
std::valarray<size_t> sort_indexes(const std::valarray<T> &v) {

    // initialize original index locations
    std::valarray<size_t> idx(v.size());
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
Array sort_by_indexes(const Array & array, const std::valarray<size_t> & indexes){
    if (array.size() != indexes.size()){
        std::cerr << AT << " size mismatch\n";
        exit(1);
    }
    Array sorted (array.size());
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
        std::cerr << AT<< "\n";
        exit(1);
    }
    if (data.empty()) {
        (*p_log)(LOG_ERR,  AT) << "Error! Empty vector of data cannot be outputed"
                          << "\n Exiting...\n";
        std::cerr << AT<< "\n";
        exit(1);
    }
    if (data[0].size() == 0) {
        (*p_log)(LOG_ERR, AT) << "Error! Empty data array cannot be outputed"
                          << "\n Exiting...\n";
        std::cerr << AT<< "\n";
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

void remove_file_if_existis(const std::string &fname){

//    LOG_INIT_CERR();
//    Log.set_log_level(LOG_SILENT);

    try {
        if (std::experimental::filesystem::remove(fname))
//            (*p_log)(LOG_INFO) << "file " << fname << " deleted.\n";
            std::cout << "file " << fname << " deleted.\n";
        else
//            (*p_log)(LOG_INFO) << "file " << fname << " not found.\n";
            std::cout << "file " << fname << " not found.\n";
    }
    catch(const std::experimental::filesystem::filesystem_error& err) {
//        (*p_log)(LOG_ERR) << "filesystem error: " << err.what() << '\n';
        std::cerr << "filesystem error: " << err.what() << '\n';
        std::cerr <<AT<<"\n";
    }
}


/*
 * Print Array in a form to copy-past in PyBlastAfterglow
 */
static void print_xy_as_numpy(double * x_arr, double * y_arr, int nxy, int ever_i=10) {
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
static void print_xy_as_numpy(Array & x_arr, Array & y_arr, int nxy, int ever_i=10) {
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
static void print_x_as_numpy(Array & x_arr, int ever_i, std::string name="x_arr", std::string label="") {
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


struct Output{
private:
//    std::unique_ptr<logger> p_log;
public:

    explicit Output(int loglevel){
//        p_log = std::make_unique<logger>(std::cout, loglevel, "Output");
    }

    template<typename vecarr>
    void Data(vecarr data,
              std::vector<std::string> &names,
              const std::string& fpath,
              const bool show_output = true){

//    LOG_INIT_CERR();
//    Log.set_log_level(LOG_SILENT);

        // remove the old file (so to avoid data piling up)
        remove_file_if_existis(fpath);

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
    void VectorOfTablesH5(std::vector<std::vector<std::vector<double>>> data,
                          std::vector<std::string> set_names,
                          std::string fpath){

//    std::cout << "S"

        // remove the old file (so to avoid data piling up)
        remove_file_if_existis(fpath);

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



    void VectorOfVectorsH5(VecVector data,
                           std::vector<std::string> set_names,
                           std::string fpath,
                           std::unordered_map<std::string, double> attrs={}){

//    std::cout << "S"
        if (data.empty()){
            std::cerr <<" Empty data cannot be saved.\n"
                      << " Exiting...\n";
            std::cerr << AT << "\n";
            exit(1);
        }
        if (data[0].empty()){
            std::cerr <<" Empty data[0] cannot be saved\n"
                      << " Exiting...\n";
            std::cerr << AT << "\n";
            exit(1);
        }
        if (data.size() != set_names.size()){
            std::cerr <<"Datasize and name.size are not the same; Data has "<<data.size()<<" names have "<<set_names.size()<<"\n"
                      << " Exiting...\n";
            std::cerr << AT << "\n";
            exit(1);
        }

        // remove the old file (so to avoid data piling up)
        remove_file_if_existis(fpath);

        H5::H5File file(fpath, H5F_ACC_TRUNC); // "/home/m/Desktop/tryout/file.h5"

        size_t n_entries = data.size();

        for (size_t id = 0; id < n_entries; ++id){

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

    void VecVectorOfVectorsAsGroupsH5(std::vector<std::vector<std::vector<double>>> data,
                                      std::vector<std::string> group_names,
                                      std::vector<std::string> arr_names,
                                      std::string fpath,
                                      std::unordered_map<std::string, double> attrs = {},
                                      std::vector<std::unordered_map<std::string, double>> group_attrs = {}
    ){
        // remove the old file (so to avoid data piling up)
        remove_file_if_existis(fpath);

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
                    std::cerr  << " size mismatch group_attrs="<<group_attrs.size()<<" group_names="<<group_names.size()<<"\n"
                               << " Exiting...\n";
                    std::cerr << AT << "\n";
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
            std::vector<std::vector<std::vector<double>>> group_data,
            std::vector<std::string> group_names,
            std::vector<std::string> vector_names_in_each_group,
            std::vector<std::vector<double>> data,
            std::vector<std::string> vector_names,
            std::unordered_map<std::string, double> attrs = {}
    ){
        // remove the old file (so to avoid group_data piling up)
        remove_file_if_existis(fpath);

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
            std::cerr <<" Empty data cannot be saved.\n"
                      << " Exiting...\n";
            std::cerr << AT << "\n";
            exit(1);
        }
        if (data[0].empty()){
            std::cerr <<" Empty data[0] cannot be saved\n"
                      << " Exiting...\n";
            std::cerr << AT << "\n";
            exit(1);
        }
        if (data.size() != vector_names.size()){
            std::cerr <<"Datasize and name.size are not the same; Data has "
                      <<data.size()<<" names have "<<vector_names.size()<<"\n"
                      << " Exiting...\n";
            std::cerr << AT << "\n";
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
            std::vector<std::vector<std::vector<std::vector<double>>>> group_data,
            std::vector<std::string> group_names_tables,
            std::vector<std::string> table_names_in_each_group,
            std::vector<std::vector<double>> data = {},
            std::vector<std::string> vector_names = {},
            std::unordered_map<std::string, double> attrs = {}
    ){
        // remove the old file (so to avoid group_data piling up)
        remove_file_if_existis(fpath);

        H5::H5File file(fpath, H5F_ACC_TRUNC); // "/home/m/Desktop/tryout/file.h5"

        /// save group data
        if (group_names_tables.size()!=group_data.size()){
            std::cerr  << " Group names="<<group_names_tables.size()
                       <<" while data vec contains=" << group_data.size() << " groups of tables" << "\n"
                       <<" Exiting...";
            std::cerr << AT << "\n";
            exit(1);
        }
        if (table_names_in_each_group.size()!=group_data[0].size()){
            std::cerr  << " Table names="<<table_names_in_each_group.size()
                       <<" while each group contains=" << group_data[0].size() << " tables" << "\n"
                       << " Exiting...\n";
            std::cerr << AT << "\n";
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

        for (size_t id = 0; id < n_entries; ++id){
            H5::Group grp(file.createGroup(group_names_tables[id]));
            for (size_t iv = 0; iv < group_data[id].size(); iv++){
                size_t nrow = group_data[id][iv].size();
                size_t ncol = group_data[id][iv][0].size();

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
            std::cerr  << " Empty data cannot be saved.\n";
            std::cerr << AT << "\n";
            exit(1);
        }
        if (data[0].empty()) {
            std::cerr  << " Empty data[0] cannot be saved\n";
            std::cerr << AT << "\n";
            exit(1);
        }
        if (data.size() != vector_names.size()) {
            std::cerr  << "Datasize and name.size are not the same; Data has " << data.size() << " names have "
                       << vector_names.size() << "\n";
            std::cerr << AT << "\n";
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



#endif //SRC_UTILS_H
