#include <iostream>
#include <fstream>
#include <random>
#include <utility>
#include <vector>
#include <random>
#include <map>
#include <algorithm>
#include <cstring>
#include <climits>
#include <unordered_set>
#include <omp.h>
#include <sys/types.h>
#include <dirent.h>
#include <boost/filesystem.hpp>

#define COUNTER_DETERMINISTIC 0
#define COUNTER_NON_DETERMINISTIC 1
#define RANDOM_RESET 2
#define TIMESTAMP 3

#define D_SIZE 10

#define DATASET "canterbury"

#define REPETITIONS 10

// 0 stands for static list configuration
// 1 stands for dynamic list configuration
#define CONFIGURATION 0

// If value of BWT is 0 do not BWT transform the files
// Otherwise BWT transform the files
#define BWT 0

// 0 stands to run the cost experiments
// 1 stands to run the opt-differences experiment
#define EXPERIMENTS 0

// Uncomment if want to run in parallel
//#define RUN_IN_PARALLEL
//#define NUMBER_THREADS 4


// Returns a vector of chars of the request sequence at file located in path
std::vector<char> read_request_sequence(const std::string &path)
{
    //Read the request sequence from the file and save it in a vector of chars;
    char ch;
    std::vector<char> request_sequence;
    std::fstream fin(path, std::fstream::in);
    while (fin >> std::noskipws >> ch) {
        request_sequence.push_back(ch);
    }
    return request_sequence;
}

// Returns a vector with all possible permutations of the elements of vector l0
std::vector<std::vector<char>> create_permutations(std::vector<char> l0)
{
    std::vector<std::vector<char>> permutations;

    std::sort(l0.begin(), l0.end());
    do {
        permutations.push_back(l0);
    } while(std::next_permutation(l0.begin(), l0.end()));

    return permutations;
}

//The cost of converting l1 to l2
//Using this algo: https://stackoverflow.com/questions/35942512/counting-inversions-between-2-arrays
int cost_move(std::vector<char> l1, std::vector<char> l2)
{
    //Create a set mapping elements of l1 to their indices
    std::map<char, int> H;
    for (int i = 0; i < l1.size(); i++)
    {
        std::pair<char, int> pair = std::make_pair(l1[i], i+1);
        H.insert(pair);
    }

    //Convert l2 accordingly
    std::vector<int> l2_int;
    for (char i : l2) {
        l2_int.push_back(H[i]);
    }

    //Find the number of inversions of l2_int
    int inversions_count = 0;

    for (int i = 0; i < l2.size()-1; i++)
    {
        for (int j = i+1; j < l2.size(); j++)
        {
            if (l2_int[i] > l2_int[j])
                inversions_count++;
        }
    }

    return inversions_count;
}

// Returns position of elemen ch in list list
int pos(char ch, std::vector<char> list)
{
    for (int pos = 0; pos < list.size(); pos++)
    {
        if (list[pos] == ch)
            return pos+1;
    }
    std::cout<<"ERROR: char not found in list"<<std::endl;

    return -1;
}

// Returns a vector with the count of the first element of the pair 'pair'
// at every position in the request sequence.
std::vector<int> find_x_counts(std::vector<char> request_sequence, std::pair<char, char> pair)
{
    size_t n = request_sequence.size();
    std::vector<int> x_counts;

    int cur_count = 0;
    x_counts.push_back(0);
    int it = 1;

    for (char elem : request_sequence)
    {
        if (elem == pair.first)
            cur_count++;
        x_counts.push_back(cur_count);
        it++;
    }
    return x_counts;
}

// Find the position of the element ch at the pairwise list pair
int find_pos_list(std::pair<char, char> pair, char ch)
{
    if (pair.first == ch)
        return 0;
    else
        return 1;
}

// Find the position of the element ch at the inverted pairwise list pair
int find_pos_list_inv(std::pair<char, char> pair, char ch)
{
    if (pair.first == ch)
        return 1;
    else
        return 0;
}

// Returns the minimum between a and b
int min(int a, int b)
{
    if (a < b)
        return a;
    return b;
}

// Returns the minimum between a, b and c
int min3(int a, int b, int c)
{
    return min(min(a, b), c);
}

// Pairwise optimal algorithm 1 - dynamic programming approach
// Running time complexity: quadratic in length of request sequence
int calculate_opt_xy(const std::vector<char> &request_sequence, std::pair<char, char> pair, int d)
{
    size_t n = request_sequence.size();
    std::vector<int> x_counts = find_x_counts(request_sequence, pair);

    std::vector<int> c1(n, INT_MAX);
    std::vector<int> c2(n, INT_MAX);

    c1[0] = min(find_pos_list(pair, request_sequence[0]), d + find_pos_list_inv(pair, request_sequence[0]) + d);
    c2[0] = min(find_pos_list_inv(pair, request_sequence[0]) + d, find_pos_list(pair, request_sequence[0]) + d);

    for (int j = 0; j < n; j++)
    {
        int min_val_c1 = INT_MAX;
        int min_val_c2 = INT_MAX;

        for (int k = 0; k < j; k++)
        {
            int cost_yx = x_counts[j+1] - x_counts[k+1];
            int cost_xy = j - k - cost_yx;

            c1[j] = min3(min_val_c1, c1[k] + cost_xy, c2[k] + cost_yx + d);
            min_val_c1 = c1[j];
            c2[j] = min3(min_val_c2, c2[k] + cost_yx, c1[k] + cost_xy + d);
            min_val_c2 = c2[j];
        }
    }
    return min(c1[n-1], c2[n-1]);
}

// Pairwise optimal algorithm 2 - running time complexity: linear in length of request sequence
int calculate_opt_xy_2(std::vector<char> request_sequence, std::pair<char, char> pair, int d)
{
    auto n = request_sequence.size();
    int cur_list_config = 0;
    int cost = 0;

    int i = 0;

    while (i < n)
    {
        char r = request_sequence[i];
        if (cur_list_config == 0)
        {
            //Here a phase starts. we iterate till we reach the end or till we find
            // a point where the number of 'y's - number of 'x's is bigger or equal than 2*d
            if (r == pair.second)
            {
                int cur_x_count = 0;
                int cur_y_count = 0;

                int j = i;

                bool do_not_cont = false;
                while (j < n)
                {
                    r = request_sequence[j];
                    if (r == pair.first)
                        cur_x_count++;
                    else
                        cur_y_count++;
                    if (cur_x_count > cur_y_count)
                    {
                        do_not_cont = true;
                        cost += cur_y_count;
                        break;
                    }
                    if (cur_y_count - cur_x_count >= 2*d)
                        break;
                    j+=1;
                }
                if (j != n && !do_not_cont)
                {
                    cost += d + cur_x_count;
                    cur_list_config = 1;
                }
                else if (!do_not_cont)
                {
                    if (cur_y_count - cur_x_count > d)
                        cost += d + cur_x_count;
                    else
                        cost += cur_y_count;
                }
                i = j+1;
            }
            else
                i++;
        }
        else if (cur_list_config == 1)
        {
            if (r == pair.first)
            {
                int cur_x_count = 0;
                int cur_y_count = 0;

                int j = i;

                bool do_not_cont = false;

                while (j < n)
                {
                    r = request_sequence[j];
                    if (r == pair.first)
                        cur_x_count += 1;
                    else
                        cur_y_count += 1;
                    if (cur_y_count > cur_x_count)
                    {
                        do_not_cont = true;
                        cost += cur_x_count;
                        break;
                    }
                    if (cur_x_count - cur_y_count >= 2*d)
                        break;
                    j++;
                }
                if (j != n && !do_not_cont)
                {
                    cost += d + cur_y_count;
                    cur_list_config = 0;
                }
                else if (!do_not_cont)
                {
                    if (cur_x_count - cur_y_count > d)
                        cost += d + cur_y_count;
                    else
                        cost += cur_x_count;
                }
                i = j+1;
            }
            else
                i++;
        }
    }
    return cost;

}

// Returns the optimal cost to serve the request sequence sigma
// with paid exchange factor d and initial list configuration l0
int real_opt_cost(const std::vector<char> &l0, const std::vector<char> &sigma, int d)
{
    size_t m = sigma.size();

    //Create all possible permutations of list l0
    std::vector<std::vector<char>> list_of_lists = create_permutations(l0);
    std::vector<std::vector<int>> dyn;

    for (const auto &list_of_list : list_of_lists) {
        std::vector<int> cost_per_request;
        cost_per_request.push_back(d*cost_move(l0, list_of_list));
        dyn.push_back(cost_per_request);
    }

    for (int j = 1; j < sigma.size()+1; j++)
    {
#ifdef RUN_IN_PARALLEL
    #pragma omp parallel for schedule(dynamic)
#endif
        for (int i = 0; i < list_of_lists.size(); i++)
        {
            //Calculate min
            int min = INT_MAX;
            for (int k = 0; k < list_of_lists.size(); k++)
            {
                int prev_dyn_cost = dyn[k][j-1];
                int pos_cost = pos(sigma[j-1], list_of_lists[k]);
                int move_cost = d*cost_move(list_of_lists[k], list_of_lists[i]);
                int cost = prev_dyn_cost + pos_cost + move_cost;
                if (cost < min)
                    min = cost;
            }
            dyn[i].push_back(min);
        }
    }

    int min = INT_MAX;
    for (int k = 0; k < list_of_lists.size(); k++)
    {
        int cost = dyn[k][m];
        if (cost < min)
            min = cost;
    }
    return min;
}

// Calculates the pairwise approximation cost of serving a request sequence
// sigma using l0 as initial list configuration and having a paid exchange factor d
int pairwise_opt_cost(const std::vector<char> &l0, std::vector<char> sigma, int d)
{
    std::vector<std::pair<char, char>> pairs;
    for (int i = 0; i < l0.size(); i++)
    {
        for (int j = i+1; j < l0.size(); j++)
        {
            pairs.emplace_back(l0[i], l0[j]);
        }
    }

    std::unordered_set<char> sequence_set;
    for (char elem : sigma)
    {
        sequence_set.insert(elem);
    }

    //Calculate all sequences sigma_xy
    std::vector < std::pair < std::vector<char>, std::pair < char, char > > > sigma_xys;
    for (auto pair : pairs) {
        std::vector<char> sigma_xy;
        if (sequence_set.count(pair.first) != 0 || sequence_set.count(pair.second) != 0)
        {
            for (char elem : sigma)
            {
                if (elem == pair.first || elem == pair.second)
                    sigma_xy.push_back(elem);
            }
            sigma_xys.emplace_back(sigma_xy, pair);
        }
    }

    int sum_cost = 0;
#ifdef RUN_IN_PARALLEL
    #pragma omp parallel for schedule(dynamic) reduction(+:sum_cost)
#endif
    for (int i = 0; i < sigma_xys.size(); i++)
    {
        std::pair<std::vector<char>, std::pair<char, char>> sigma_xy = sigma_xys[i];
        if (!sigma_xy.first.empty())
        {
            // Optionally calculate_opt_xy(sigma_xy.first, sigma_xy.second, d) can be used
            // to calculate the pairwise optimal cost using the dynamic programming approach
            int pairwise_opt = calculate_opt_xy_2(sigma_xy.first, sigma_xy.second, d);
            sum_cost += pairwise_opt;
        }
    }
    return static_cast<int>(sum_cost + sigma.size());

}

/*
 * This function reads the input sequences from the files located in dataset_path
 * and calculates the differences between the cost of the pairwise approximation
 * to serve the sequence and the cost of the real opt. It picks uniformly at random
 * 3, 4 and 5 elements from the input sequence and then extracts the sequence
 * consisting only of these elements from the original input sequence. It serves the
 * new sequences on these lists using both algorithms and writes the results in
 * a file in results_path.
 */
void calculate_opts_differences(const std::string &dataset_path, const std::string &results_path,
                    std::vector<std::string> filenames)
{
    std::ofstream results_file;

    for(const std::string &filename : filenames)
    {
        results_file.open(results_path, std::ios::app);
        results_file<<"FILE: "<<filename<<"\n";
        results_file.close();

        std::string path = dataset_path + filename;
        std::vector<char> request_sequence = read_request_sequence(path);

        //Build a hashmap of the elements in request_sequence;
        std::unordered_set<char> sequence_set;
        for (char elem : request_sequence)
        {
            sequence_set.insert(elem);
        }
        //Convert set to vector
        std::vector<char> elements_vec;
        elements_vec.insert(elements_vec.end(), sequence_set.begin(), sequence_set.end());

        int ds[D_SIZE] = {1, 2, 3, 4, 5, 6, 10, 20, 50, 100};

        //Pick at random the elements to include in the lists
        int list_sizes[] = {3, 4, 5};

        for (int list_size : list_sizes)
        {
            if (elements_vec.size() < list_size)
            {
                results_file<<"The doc contains only: "<<elements_vec.size()<<" unique elements"<<std::endl;
                break;
            }

            //Repeat REPEAT times in order for scientific correctness
            for (int repeat_ind = 0; repeat_ind < REPETITIONS; repeat_ind++)
            {
                std::vector<char> list;
                std::unordered_set<char> list_set;

                //Shuffle the elements of the vector containing the unique elements of the sequence
                std::vector<unsigned int> indices(elements_vec.size());
                std::iota(indices.begin(), indices.end(), 0);
                std::shuffle(indices.begin(), indices.end(), std::mt19937(std::random_device()()));

                for (int i = 0; i < list_size; i++)
                {
                    list.push_back(elements_vec[indices[i]]);
                    list_set.insert(elements_vec[indices[i]]);
                }


                results_file.open(results_path, std::ios::app);
                results_file<<"List size: "<<list_size<<"\nList:\t";
                for (char ch : list)
                {
                    results_file<<(int)ch<<"\t";
                }
                results_file<<"\n";
                results_file.close();

                //Calculate the sequence sigma only with the elements found above
                std::vector<char> sigma;
                for (char ch : request_sequence)
                {
                    if (list_set.count(ch) != 0)
                        sigma.push_back(ch);
                }

                results_file.open(results_path, std::ios::app);
                results_file<<"Sigma length: "<<sigma.size()<<"\n";
                results_file.close()    ;

                results_file.open(results_path, std::ios::app);
                for (int d : ds)
                {
                    //Calculate cost here
                    int sum_real_cost = real_opt_cost(list, sigma, d);
                    int sum_pairwise_cost = pairwise_opt_cost(list, sigma, d);

                    //Write cost in the file
                    results_file<<d<<"\t"<<sum_real_cost<<"\t"<<sum_pairwise_cost<<"\t"<<sum_pairwise_cost - sum_real_cost<<std::endl;
                }
                results_file<<"\n";
                results_file.close();
            }

        }
    }
}

// Generates the list consisting of all the ASCII characters with
// counters initialized to 0
std::vector<std::pair<char, int>> generate_counter_deterministic_list()
{
    std::vector<std::pair<char, int>> list;

    //Start from character '0' and go till the 256th character in ASCII
    char ch='0';
    int cnt =0;
    while (cnt < 256)
    {
        list.emplace_back(ch, 0);
        cnt++;
        ch++;
    }
    return list;
}

// Generates the list consisting of all the ASCII characters with
// counters initilized uniformly at random with a number 0 to k-1
std::vector<std::pair<char, int>> generate_list_counter_uniformly(int k)
{
    std::vector<std::pair<char, int>> list;

    char ch='0';
    int cnt =0;

    const int range_from  = 0;
    const int range_to    = k-1;
    std::random_device                  rand_dev;
    std::mt19937                        generator(rand_dev());
    std::uniform_int_distribution<int>  distr(range_from, range_to);

    while (cnt < 256)
    {
        //Calculate counter - random number between 0 and k-1.
        int counter = distr(generator);
        list.emplace_back(ch, counter);
        cnt++;
        ch++;
    }

    return list;
}

// Generates the list consisting of all the ASCII characters with
// counters initialized according to a distribution given as input
std::vector<std::pair<char, int>> generate_list_from_distribution(std::vector<double> distribution)
{
    std::vector<std::pair<char, int>> list;

    char ch='0';
    int cnt =0;

    std::random_device                  rand_dev;
    std::mt19937                        generator(rand_dev());
    std::discrete_distribution<> distr(distribution.begin(), distribution.end());

    while (cnt < 256)
    {
        //Calculate counter - random number between 1 and k.
        int counter = distr(generator);
        list.emplace_back(ch, counter+1);
        cnt++;
        ch++;
    }

    return list;
}

/*
 * Serves the request sequence using the counter(k, {k-1}) algorithm with initial list configuration list
 * and paid exchange factor d. Returns the cost to serve the sequence in the variable total_paid_exchange_cost
 */
int counter(std::vector<char> request_sequence, std::vector<std::pair<char, int>> list, int k,
            int d, int &total_paid_exchange_cost)
{
    int total_cost = 0;
    //int total_paid_exchange_cost = 0;

    //Go through each character of the sequence
    for (char request_element : request_sequence)
    {
        int service_cost = 0;

        bool found = false;
        //Go through the list and find the element corresponding to char
        //for (std::pair<char, int> list_element : list)
        int list_index;
        for (list_index = 0; list_index < list.size(); list_index++)
        {
            service_cost++;
            if (request_element == list[list_index].first)
            {
                total_cost += service_cost;
                found = true;

                //Decrement counter - careful when counter is 0
                if (list[list_index].second == 0)
                    list[list_index].second = k-1;
                else
                    list[list_index].second--;

                //Move item to front if the counter is k-1
                if (list[list_index].second == k-1)
                {
                    //Calculate the number of paid exchanges
                    total_paid_exchange_cost += list_index*d;

                    //Shift all elements before ch one to the left and move element ch to the front of the list
                    std::pair<char, int> current_item = list[list_index];
                    for(int i = list_index; i > 0; i--)
                    {
                        list[i] = list[i-1];
                    }
                    list[0] = current_item;
                }

                break;
            }
        }
        if (!found)
        {
            // If dynamic list configuration and item not found, insert it in the list
            if (CONFIGURATION == 1)
            {
                //std::cout<<"Char "<<request_element<<" not found"<<std::endl;
                service_cost++;
                total_cost += service_cost;
                list.emplace_back(request_element, 0);

                //Decrement the counter
                if (list[list_index].second == 0)
                    list[list_index].second = k-1;
                else
                    list[list_index].second--;

                //Move item to front if the counter is k-1
                if (list[list_index].second == k-1)
                {
                    //Calculate the number of paid exchanges
                    total_paid_exchange_cost += list_index*d;

                    //Move To Front
                    //Shift all elements before ch one to the left and move element ch to the front of the list
                    std::pair<char, int> current_item = list[list_index];
                    for(int i = list_index; i > 0; i--)
                    {
                        list[i] = list[i-1];
                    }
                    list[0] = current_item;
                }
            }
            else
                std::cout<<"Char "<<request_element<<" not found"<<std::endl;
        }
    }
    return total_cost;
}

/*
 * Serves the request sequence using the random_reset(k) algorithm with initial list configuration list,
 * paid exchange factor d, a given stationary distribution and a given resetting distribution.
 * Returns the cost to serve the sequence in the variable total_paid_exchange_cost.
 */
int random_reset(std::vector<char> request_sequence, std::vector<std::pair<char, int>> list, double reset_to_k_prob,
                 std::vector<double> stationary_distribution, int k, int d, int &total_paid_exchange_cost)
{
    int total_cost = 0;

    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    std::discrete_distribution<> distr(stationary_distribution.begin(), stationary_distribution.end());

    //Go through each character of the sequence
    for (char request_element : request_sequence)
    {
        int service_cost = 0;

        bool found = false;
        //Go through the list and find the element corresponding to char
        //for (std::pair<char, int> list_element : list)
        int list_index;
        for (list_index = 0; list_index < list.size(); list_index++)
        {
            service_cost++;
            if (request_element == list[list_index].first)
            {
                total_cost += service_cost;
                found = true;

                //Decrement counter - careful when counter is 0
                if (list[list_index].second > 1)
                    list[list_index].second--;

                //Move item to front if the counter is k-1
                if (list[list_index].second == 1)
                {
                    //Calculate the number of paid exchanges
                    total_paid_exchange_cost += list_index*d;

                    //Move To Front
                    //Shift all elements before ch one to the left and move element ch to the front of the list
                    std::pair<char, int> current_item = list[list_index];
                    for(int i = list_index; i > 0; i--)
                    {
                        list[i] = list[i-1];
                    }
                    list[0] = current_item;

                    //Reset the value of the counter to k with probability reset_to_k_prob,
                    //and reset to k-1 with probability 1-reset_to_k_prob

                    std::random_device rd;
                    std::mt19937 gen(rd());
                    std::uniform_real_distribution<> distr(0.0, 1.0);
                    double result = distr(gen);
                    if(result < reset_to_k_prob)
                        list[0].second = k;
                    else
                        list[0].second = k-1;
                }

                break;
            }
        }
        if (!found)
        {

            // If dynamic list configuration - insert item in the list
            if (CONFIGURATION == 1)
            {
                int counter = distr(generator);
                list.emplace_back(request_element, counter+1);
                service_cost++;
                total_cost += service_cost;

                //Decrement counter - careful when counter is 0
                if (list[list_index].second > 1)
                    list[list_index].second--;

                //Move item to front if the counter is k-1
                if (list[list_index].second == 1)
                {
                    //Calculate the number of paid exchanges
                    total_paid_exchange_cost += list_index*d;

                    //Move To Front
                    //Shift all elements before ch one to the left and move element ch to the front of the list
                    std::pair<char, int> current_item = list[list_index];
                    for(int i = list_index; i > 0; i--)
                    {
                        list[i] = list[i-1];
                    }
                    list[0] = current_item;

                    //Reset the value of the counter to k with probability reset_to_k_prob,
                    //and reset to k-1 with probability 1-reset_to_k_prob

                    std::random_device rd;
                    std::mt19937 gen(rd());
                    std::uniform_real_distribution<> distr(0.0, 1.0);
                    double result = distr(gen);
                    if(result < reset_to_k_prob)
                        list[0].second = k;
                    else
                        list[0].second = k-1;
                }
            }
            else
                std::cout<<"Char "<<request_element<<" not found"<<std::endl;
        }
    }
    return total_cost;
}

/*
 * Serves the request sequence using the timestamp(k, p) algorithm with initial list configuration list,
 * paid exchange factor d and policy probability p.
 * Returns the cost to serve the sequence in the variable total_paid_exchange_cost.
 */
int timestamp(std::vector<char> request_sequence, std::vector<std::pair<char, int>> list, double p,
              int k, int d, int &total_paid_exchange_cost)
{
    int total_cost = 0;
    const int range_from  = 0;
    const int range_to    = k-1;
    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    std::uniform_int_distribution<int>  distr(range_from, range_to);

    //Bookkeeping for the request sequence
    //Keep a pair for each request item - a bool saying whether the master request was performed
    //and a bool saying whether Policy 2 was performed.
    std::vector<std::pair<bool, bool>> request_seq_master_policy;

    //Go through each character of the sequence
    //for (char request_element : request_sequence)
    for (int request_element_index = 0; request_element_index < request_sequence.size(); request_element_index++)
    {
        char request_element = request_sequence[request_element_index];
        int service_cost = 0;

        bool found = false;
        //Go through the list and find the element corresponding to char
        //for (std::pair<char, int> list_element : list)
        int list_index;
        for (list_index = 0; list_index < list.size(); list_index++)
        {
            service_cost++;
            if (request_element == list[list_index].first)
            {
                total_cost += service_cost;
                found = true;

                //Decrement counter - careful when counter is 0
                if (list[list_index].second == 0)
                    list[list_index].second = k-1;
                else
                    list[list_index].second--;

                //Move item to front if the counter is k-1
                if (list[list_index].second == k-1)
                {
                    std::random_device rd;
                    std::mt19937 gen(rd());
                    std::uniform_real_distribution<> distr(0.0, 1.0);
                    double result = distr(gen);

                    //With probability p move item to the front == Perform Policy 1
                    if(result < p)
                    {
                        //Calculate the number of paid exchanges
                        total_paid_exchange_cost += list_index*d;

                        //Shift all elements before ch one to the left and move element ch to the front of the list
                        std::pair<char, int> current_item = list[list_index];
                        for(int i = list_index; i > 0; i--)
                        {
                            list[i] = list[i-1];
                        }
                        list[0] = current_item;

                        //The master request to x occurs but policy 1 is performed (policy 2 is not performed)
                        request_seq_master_policy.emplace_back(true, false);
                    }

                    //With probability 1-p perform Policy 2
                    else
                    {
                        //Both master request occurs and policy 2 is performed
                        request_seq_master_policy.emplace_back(true, true);

                        //Find the longest suffix lambda(t) that ends with request_element and contains
                        //exactly one master request to request_element.

                        //Keep a dict where you save each element traversed while finding lambda and
                        //the number of master requests to this element in lambda.
                        std::map<char, std::pair<int, bool>> element_master_requests_map;

                        int suffix_index;
                        for (suffix_index = request_element_index - 1; suffix_index >= 0; suffix_index--)
                        {
                            //Here you find the first element that is equal to x and it has been a master request
                            if(request_sequence[suffix_index] == request_element
                               && request_seq_master_policy[suffix_index].first == 1)
                            {
                                suffix_index++;
                                break;
                            }
                            else
                            {
                                if (element_master_requests_map.count(request_sequence[suffix_index]) == 0)
                                {
                                    element_master_requests_map.insert(std::make_pair(request_sequence[suffix_index],
                                                                                      std::make_pair(request_seq_master_policy[suffix_index].first,
                                                                                                     request_seq_master_policy[suffix_index].second)));
                                }
                                else
                                {
                                    element_master_requests_map[request_sequence[suffix_index]].first +=
                                            request_seq_master_policy[suffix_index].first;
                                    element_master_requests_map[request_sequence[suffix_index]].second =
                                            element_master_requests_map[request_sequence[suffix_index]].second
                                            || request_seq_master_policy[suffix_index].second;

                                }
                            }
                        }
                        //suffix_index now is the index of the first element of the sequence lambda(t)
                        //request_element_index is the index of the last element of the sequence lambda(t)


                        //Find the first item in the list (z) such that there is at most one master request to z
                        //in lambda(t) and the request was served using Policy 2.
                        //Go through all elements in the list
                        int move_to_index;
                        for (int list_index_2 = 0; list_index_2 < list.size(); list_index_2++)
                        {
                            char list_element = list[list_index_2].first;

                            //Check if there is at most one master request to list_element in lambda(t)
                            if (element_master_requests_map.count(list_element))
                            {
                                if (element_master_requests_map[list_element].first == 1 &&
                                    element_master_requests_map[list_element].second)
                                {
                                    move_to_index = list_index_2;
                                    break;
                                }
                                else if (element_master_requests_map[list_element].first == 0)
                                {
                                    move_to_index = list_index_2;
                                    break;
                                }
                            }
                            else
                            {
                                move_to_index = list_index_2;
                                break;
                            }
                        }

                        //Move the current element before the element list[move_to_index]
                        //Increase paid exchange cost
                        std::pair<char, int> current_item = list[list_index];
                        for(int i = list_index; i > move_to_index; i--)
                        {
                            total_paid_exchange_cost+=d;
                            list[i] = list[i-1];
                        }
                        total_paid_exchange_cost+=d;
                        list[move_to_index] = current_item;
                    }
                }

                else
                {
                    //Neither master request occurs, nor policy 2 is performed (policy 1 is performed)
                    request_seq_master_policy.emplace_back(false, false);
                }

                break;
            }
        }
        if (!found)
        {
            // If dynamic list configuration
            if (CONFIGURATION == 1)
            {
                int counter = distr(generator);
                list.emplace_back(request_element, counter);

                service_cost++;
                total_cost += service_cost;


                if (list[list_index].second == 0)
                    list[list_index].second = k-1;
                else
                    list[list_index].second--;

                //Move item to front if the counter is k-1
                if (list[list_index].second == k-1)
                {
                    std::random_device rd;
                    std::mt19937 gen(rd());
                    std::uniform_real_distribution<> distr(0.0, 1.0);
                    double result = distr(gen);

                    //With probability p move item to the front == Perform Policy 1
                    if(result < p)
                    {
                        //Calculate the number of paid exchanges
                        total_paid_exchange_cost += list_index*d;

                        //Shift all elements before ch one to the left and move element ch to the front of the list
                        std::pair<char, int> current_item = list[list_index];
                        for(int i = list_index; i > 0; i--)
                        {
                            list[i] = list[i-1];
                        }
                        list[0] = current_item;

                        //The master request to x occurs but policy 1 is performed (policy 2 is not performed)
                        request_seq_master_policy.emplace_back(true, false);
                    }

                        //With probability 1-p perform Policy 2
                    else
                    {
                        //Both master request occurs and policy 2 is performed
                        request_seq_master_policy.emplace_back(true, true);

                        //Find the longest suffix lambda(t) that ends with request_element and contains
                        //exactly one master request to request_element.

                        //Keep a dict where you save each element traversed while finding lambda and
                        //the number of master requests to this element in lambda.
                        std::map<char, std::pair<int, bool>> element_master_requests_map;

                        int suffix_index;
                        for (suffix_index = request_element_index - 1; suffix_index >= 0; suffix_index--)
                        {
                            //Here you find the first element that is equal to x and it has been a master request
                            if(request_sequence[suffix_index] == request_element
                               && request_seq_master_policy[suffix_index].first == 1)
                            {
                                suffix_index++;
                                break;
                            }
                            else
                            {
                                if (element_master_requests_map.count(request_sequence[suffix_index]) == 0)
                                {
                                    element_master_requests_map.insert(std::make_pair(request_sequence[suffix_index],
                                                                                      std::make_pair(request_seq_master_policy[suffix_index].first,
                                                                                                     request_seq_master_policy[suffix_index].second)));
                                }
                                else
                                {
                                    element_master_requests_map[request_sequence[suffix_index]].first +=
                                            request_seq_master_policy[suffix_index].first;
                                    element_master_requests_map[request_sequence[suffix_index]].second =
                                            element_master_requests_map[request_sequence[suffix_index]].second
                                            || request_seq_master_policy[suffix_index].second;

                                }
                            }
                        }
                        //suffix_index now is the index of the first element of the sequence lambda(t)
                        //request_element_index is the index of the last element of the sequence lambda(t)


                        //Find the first item in the list (z) such that there is at most one master request to z
                        //in lambda(t) and the request was served using Policy 2.
                        //Go through all elements in the list
                        int move_to_index;
                        for (int list_index_2 = 0; list_index_2 < list.size(); list_index_2++)
                        {
                            char list_element = list[list_index_2].first;

                            //Check if there is at most one master request to list_element in lambda(t)
                            if (element_master_requests_map.count(list_element))
                            {
                                if (element_master_requests_map[list_element].first == 1 &&
                                    element_master_requests_map[list_element].second)
                                {
                                    move_to_index = list_index_2;
                                    break;
                                }
                                else if (element_master_requests_map[list_element].first == 0)
                                {
                                    move_to_index = list_index_2;
                                    break;
                                }
                            }
                            else
                            {
                                move_to_index = list_index_2;
                                break;
                            }
                        }

                        //Move the current element before the element list[move_to_index]
                        //Increase paid exchange cost
                        std::pair<char, int> current_item = list[list_index];
                        for(int i = list_index; i > move_to_index; i--)
                        {
                            total_paid_exchange_cost+=d;
                            list[i] = list[i-1];
                        }
                        total_paid_exchange_cost+=d;
                        list[move_to_index] = current_item;
                    }
                }

                else
                {
                    //Neither master request occurs, nor policy 2 is performed (policy 1 is performed)
                    request_seq_master_policy.emplace_back(false, false);
                }
            }
            else
                 std::cout<<"Char "<<request_element<<" not found"<<std::endl;
        }
    }
    return total_cost;
}


/*
 * This function serves the input request sequence given as input with the algorithm given as input.
 * It writes the results in a file in results_path.
 */
void serve_request(const std::vector<char> &request_sequence, int algorithm, std::vector<std::pair<int, int>> d_and_ks,
                   const std::vector<double> &reset_distribution,
                   const std::vector<std::vector<double>> &stationary_distributions,
                   const std::string &results_path)
{
    std::vector<std::pair<int, int>> total_cost;
    std::ofstream results_file;
    results_file.open(results_path, std::ios::app);

    if (algorithm == COUNTER_DETERMINISTIC)
    {
        results_file<<"ALGORITHM: COUNTER DETERMINISTIC"<<std::endl;
        for (std::pair<int, int> d_and_k : d_and_ks)
        {
            int d = d_and_k.first;
            int k = d_and_k.second;
            int paid_exchange_cost = 0;
            std::vector<std::pair<char, int>> counter_list_deterministic;
            if (CONFIGURATION == 0)
                counter_list_deterministic = generate_counter_deterministic_list();
            int service_cost = counter(request_sequence, counter_list_deterministic, k, d, paid_exchange_cost);
            total_cost.emplace_back(service_cost, paid_exchange_cost);
        }
    }

    else if (algorithm == COUNTER_NON_DETERMINISTIC)
    {
        results_file<<"ALGORITHM: COUNTER NON DETERMINISTIC"<<std::endl;
        for (std::pair<int, int> d_and_k : d_and_ks)
        {
            int d = d_and_k.first;
            int k = d_and_k.second;
            std::vector<std::pair<char, int>> counter_list_uniformly;
            if (CONFIGURATION == 0)
                counter_list_uniformly = generate_list_counter_uniformly(k);
            int paid_exchange_cost = 0;
            int service_cost = counter(request_sequence, counter_list_uniformly, k, d, paid_exchange_cost);
            total_cost.emplace_back(service_cost, paid_exchange_cost);
        }
    }

    else if (algorithm == RANDOM_RESET)
    {
        results_file<<"ALGORITHM: RANDOM RESET"<<std::endl;
        int cnt =0;
        for (std::pair<int, int> d_and_k : d_and_ks)
        {
            int d = d_and_k.first;
            int k = d_and_k.second;

            std::vector<std::pair<char, int>> initial_list;

            if (CONFIGURATION == 0)
                initial_list = generate_list_from_distribution(stationary_distributions[cnt]);

            int paid_exchange_cost = 0;
            int service_cost = random_reset(request_sequence, initial_list,
                                            reset_distribution[cnt], stationary_distributions[cnt],
                                            k, d, paid_exchange_cost);
            total_cost.emplace_back(service_cost, paid_exchange_cost);
            cnt++;
        }
    }

    else if (algorithm == TIMESTAMP)
    {
        results_file<<"ALGORITHM: TIMESTAMP"<<std::endl;
        int cnt = 0;
        for (std::pair<int, int> d_and_k : d_and_ks)
        {
            int d = d_and_k.first;
            int k = d_and_k.second;
            std::vector<std::pair<char, int>> list;

            if (CONFIGURATION == 0)
                list = generate_list_counter_uniformly(k);

            int paid_exchange_cost = 0;
            int service_cost = timestamp(request_sequence, list, reset_distribution[cnt],
                                         k, d, paid_exchange_cost);
            total_cost.emplace_back(service_cost, paid_exchange_cost);
            cnt++;
        }
    }

    //Print for all different values of k and d the service cost and the paid exchange cost

    results_file<<"d\tk\tService cost\tPaid exchange cost"<<std::endl;
    for (int i = 0; i < d_and_ks.size(); i++)
    {
        results_file<<d_and_ks[i].first<<"\t";
        results_file<<d_and_ks[i].second<<"\t";
        results_file<<total_cost[i].first<<"\t";
        results_file<<total_cost[i].second<<"\n";
    }
    results_file<<"\n\n";
    results_file.close();

}

// Utility struct to build the suffix array
// The code to build the suffix array is taken from: https://codeforces.com/blog/entry/4025
struct Comparator {
    Comparator(int* pos, int gap, int N) : pos(pos), gap(gap), N(N) {};

    int* pos;
    int gap;
    int N;

    bool operator()(int i, int j) {
        if (pos[i] != pos[j])
            return pos[i] < pos[j];
        i += gap;
        j += gap;
        bool return_value = (i < N && j < N) ? pos[i] < pos[j] : i > j;
        return return_value;
    }
};

// Builds the suffix array of a requence sequence provided as input
void buildSuffixArray(std::vector<char> input, int* suffix_array)
{
    size_t N = input.size();
    auto * pos = new int[N];
    auto * tmp = new int[N];

    for (int i =0; i < N; i++)
    {
        suffix_array[i] = i;
        pos[i] = input[i];
        tmp[i] = 0;
    }
    for (int gap = 1;; gap *= 2)
    {
        std::sort(suffix_array, suffix_array + N, Comparator(pos, gap, N));
        for (int i = 0; i < (int)(N - 1); ++i)
            tmp[i + 1] = tmp[i] + Comparator(pos, gap, N)(suffix_array[i], suffix_array[i + 1]);
        for (int i = 0; i < (int)(N); ++i)
            pos[suffix_array[i]] = tmp[i];
        if (tmp[N - 1] == N - 1) break;
    }
}

// Computes the BWT array from the suffix array and the input
std::vector<char> compute_bwt_array(std::vector<int> suffix_array, std::vector<char> input)
{
    size_t n = suffix_array.size();

    if (suffix_array.size() != input.size()/2)
        std::cout<<"The sizes of suffix array and input are different"<<std::endl;

    std::vector<char> bwt_array;

    for (int i = 0; i < n; i++)
    {
        bwt_array.push_back(input[(suffix_array[i] - 1 + n) % n]);
    }
    return bwt_array;
}

// Returns the BWT transformed sequence of the input
std::vector<char> bwt(const std::vector<char> &input)
{
    // Builds the suffix array of the input
    auto* suffix_array = new int[input.size()];
    buildSuffixArray(input, suffix_array);

    size_t n = input.size()/2;
    std::vector<int> sa(n);
    int sa_index = 0;
    for (int i = 0; i < input.size(); i++)
    {
        if (suffix_array[i] < n)
        {
            sa[sa_index] = suffix_array[i];
            sa_index++;
        }
    }

    return compute_bwt_array(sa, input);
}

/*
 * This is a utility function to run all algorithms on the input sequences located in the
 * files of the folder dataset_path. It reads the input sequences, and runs all algorithms
 * and calculates the cost of serving each of these input sequences.
 * The algorithms are run 10 (REPETITIONS) times for scientific correctness
 * It writes the results in a file located in results_path
 */
void run_all_algorithms(const std::string &dataset_path, std::string& results_path,
                        std::vector<std::string> filenames)
{
    for (int iter = 1; iter <= REPETITIONS; iter++)
    {
        for(const std::string &filename : filenames)
        {
            std::string path = dataset_path + filename;
            results_path = results_path + std::to_string(iter) + ".txt";

            std::ofstream results_file;
            results_file.open(results_path, std::ios::app);
            results_file<<"FILE: "<<filename<<"\n";
            results_file.close();


            std::vector<char> request_sequence = read_request_sequence(path);

            if (BWT)
            {
                request_sequence.insert(request_sequence.end(),
                                        request_sequence.begin(),
                                        request_sequence.end());
                request_sequence = bwt(request_sequence);
            }

            std::vector<double> prob_values;
            std::vector<std::vector<double>> prob_values_2;
            //COUNTER DETERMINISTIC
            //Generate the list of d and respective k for deterministic algorithm
            std::vector<std::pair<int, int>> d_and_k_deterministic;
            d_and_k_deterministic.emplace_back(1, 1);
            d_and_k_deterministic.emplace_back(2, 2);
            d_and_k_deterministic.emplace_back(3, 2);
            d_and_k_deterministic.emplace_back(4, 3);
            d_and_k_deterministic.emplace_back(5, 4);
            d_and_k_deterministic.emplace_back(6, 5);
            d_and_k_deterministic.emplace_back(10, 8);
            d_and_k_deterministic.emplace_back(20, 16);
            d_and_k_deterministic.emplace_back(50, 39);
            d_and_k_deterministic.emplace_back(100, 78);

            serve_request(request_sequence, 0, d_and_k_deterministic, prob_values,
                          prob_values_2, results_path);


            //COUNTER NON DETERMINISTIC
            std::vector<std::pair<int, int>> d_and_k_non_deterministic;
            d_and_k_non_deterministic.emplace_back(1, 2);
            d_and_k_non_deterministic.emplace_back(2, 5);
            d_and_k_non_deterministic.emplace_back(3, 7);
            d_and_k_non_deterministic.emplace_back(4, 10);
            d_and_k_non_deterministic.emplace_back(5, 12);
            d_and_k_non_deterministic.emplace_back(6, 15);
            d_and_k_non_deterministic.emplace_back(10, 25);
            d_and_k_non_deterministic.emplace_back(20, 51);
            d_and_k_non_deterministic.emplace_back(50, 128);
            d_and_k_non_deterministic.emplace_back(100, 256);

            serve_request(request_sequence, 1, d_and_k_non_deterministic, prob_values, prob_values_2,
                          results_path);

            //RANDOM RESET
            std::vector<std::pair<int, int>> d_and_k_random_reset;
            d_and_k_random_reset.emplace_back(1, 3);
            d_and_k_random_reset.emplace_back(2, 5);
            d_and_k_random_reset.emplace_back(3, 8);
            d_and_k_random_reset.emplace_back(4, 10);
            d_and_k_random_reset.emplace_back(5, 13);
            d_and_k_random_reset.emplace_back(6, 15);
            d_and_k_random_reset.emplace_back(10, 26);
            d_and_k_random_reset.emplace_back(20, 51);
            d_and_k_random_reset.emplace_back(50, 128);
            d_and_k_random_reset.emplace_back(100, 256);

            std::vector<double> prob_values_random_reset;
            prob_values_random_reset.push_back(0.215);
            prob_values_random_reset.push_back(0.760);
            prob_values_random_reset.push_back(0.314);
            prob_values_random_reset.push_back(0.878);
            prob_values_random_reset.push_back(0.433);
            prob_values_random_reset.push_back(1);
            prob_values_random_reset.push_back(0.240);
            prob_values_random_reset.push_back(0.854);
            prob_values_random_reset.push_back(0.699);
            prob_values_random_reset.push_back(0.777);


            std::vector<double> second_to_last_prob;
            second_to_last_prob.push_back(0.451);
            second_to_last_prob.push_back(0.210);
            second_to_last_prob.push_back(0.137);
            second_to_last_prob.push_back(0.101);
            second_to_last_prob.push_back(0.080);
            second_to_last_prob.push_back(1.0/15);
            second_to_last_prob.push_back(0.04);
            second_to_last_prob.push_back(0.02);
            second_to_last_prob.push_back(0.008);
            second_to_last_prob.push_back(0.004);

            std::vector<double> last_prob;
            last_prob.push_back(0.0971);
            last_prob.push_back(0.16);
            last_prob.push_back(0.043);
            last_prob.push_back(0.089);
            last_prob.push_back(0.35);
            last_prob.push_back(1.0/15);
            last_prob.push_back(0.0095);
            last_prob.push_back(0.017);
            last_prob.push_back(0.005);
            last_prob.push_back(0.003);

            std::vector<std::vector<double>> stationary_distributions;
            for (std::pair<int, int> d_and_k : d_and_k_random_reset)
            {
                std::vector<double> stationary_distribution;
                int k = d_and_k.second;
                int cnt = 0;
                for (int j = 0; j < k-1; j++)
                {
                    stationary_distribution.push_back(second_to_last_prob[cnt]);
                }
                stationary_distribution.push_back(last_prob[cnt]);

                stationary_distributions.push_back(stationary_distribution);
            }

            serve_request(request_sequence, 2, d_and_k_random_reset, prob_values_random_reset,
                          stationary_distributions, results_path);


            //TIMESTAMP
            std::vector<std::pair<int, int>> d_and_k_timestamp;
            d_and_k_timestamp.emplace_back(1, 1);
            d_and_k_timestamp.emplace_back(2, 2);
            d_and_k_timestamp.emplace_back(3, 4);
            d_and_k_timestamp.emplace_back(4, 5);
            d_and_k_timestamp.emplace_back(5, 6);
            d_and_k_timestamp.emplace_back(6, 7);
            d_and_k_timestamp.emplace_back(10, 12);
            d_and_k_timestamp.emplace_back(20, 26);
            d_and_k_timestamp.emplace_back(50, 60);
            d_and_k_timestamp.emplace_back(100, 119);

            std::vector<double> prob_values_timestamp;
            prob_values_timestamp.push_back(0.45793);
            prob_values_timestamp.push_back(0.45798);
            prob_values_timestamp.push_back(0.4582);
            prob_values_timestamp.push_back(0.4579);
            prob_values_timestamp.push_back(0.458);
            prob_values_timestamp.push_back(0.45787);
            prob_values_timestamp.push_back(0.38835);
            prob_values_timestamp.push_back(0.5);
            prob_values_timestamp.push_back(0.459);
            prob_values_timestamp.push_back(0.4579);


            serve_request(request_sequence, 3, d_and_k_timestamp, prob_values_timestamp,
                          prob_values_2, results_path);
        }
    }
}

// Utility function to create a directory in folder_path
void create_directory_boost(const std::string &folder_path)
{
    const char* path = folder_path.c_str();
    boost::filesystem::path dir(path);
    boost::filesystem::create_directory(dir);
}

void create_directories()
{
    //Create results folders
    std::string results_path_folder = "./../results/";
    create_directory_boost("./../results/");
    create_directory_boost(results_path_folder + "opt_differences/");

    std::string results_costs_path_folder = results_path_folder + "costs/";
    create_directory_boost(results_costs_path_folder);

    std::string resutls_costs_static_path_folder = results_costs_path_folder + "static/";
    create_directory_boost(resutls_costs_static_path_folder);
    create_directory_boost(resutls_costs_static_path_folder + "canterbury/");
    create_directory_boost(resutls_costs_static_path_folder + "calgary/");
    create_directory_boost(resutls_costs_static_path_folder + "large/");

    std::string resutls_costs_dynamic_path_folder = results_costs_path_folder + "dynamic/";
    create_directory_boost(resutls_costs_dynamic_path_folder);
    create_directory_boost(resutls_costs_dynamic_path_folder + "canterbury/");
    create_directory_boost(resutls_costs_dynamic_path_folder + "calgary/");
    create_directory_boost(resutls_costs_dynamic_path_folder + "large/");

    std::string resutls_costs_bwt_path_folder = results_costs_path_folder + "bwt/";
    create_directory_boost(resutls_costs_bwt_path_folder);
    create_directory_boost(resutls_costs_bwt_path_folder + "canterbury/");
    create_directory_boost(resutls_costs_bwt_path_folder + "calgary/");
    create_directory_boost(resutls_costs_bwt_path_folder + "large/");
}

int main() {

    // There needs to be a folder named datasets inside which, there needs to be a folder
    // for each of the corpora with the respective files inside

    //The path where the folder of the datasets is located
    std::string dataset_path_folder = "./../datasets/";
    std::string dataset_path = dataset_path_folder + DATASET + "/";

    //Contains all filenames of the current dataset in use
    std::vector<std::string> filenames;

    //Read all the files in the directory of the dataset
    DIR* dirp = opendir(dataset_path.c_str());
    struct dirent * dp;
    while ((dp = readdir(dirp)) != nullptr) {
        filenames.emplace_back(dp->d_name);
    }
    closedir(dirp);

    //Create all the directories where the results will be stored
    create_directories();

#ifdef RUN_IN_PARALLEL
    omp_set_num_threads(NUMBER_THREADS);
#endif

    // If 0, run the costs experiments
    // otherwise run the opt-differences experiments
    if (EXPERIMENTS == 0)
    {
        if (BWT)
        {
            std::string results_path_folder = "./../results/costs/bwt/";
            std::string results_path = results_path_folder + DATASET + "/" + DATASET;
            run_all_algorithms(dataset_path, results_path, filenames);
        }
        else if (CONFIGURATION == 0)
        {
            std::string results_path_folder = "./../results/costs/static/";
            std::string results_path = results_path_folder + DATASET + "/" + DATASET;
            run_all_algorithms(dataset_path, results_path, filenames);
        }
        else if (CONFIGURATION == 1)
        {
            std::string results_path_folder = "./../results/costs/dynamic/";
            std::string results_path = results_path_folder + DATASET + "/" + DATASET;
            run_all_algorithms(dataset_path, results_path, filenames);
        }
    }
    else if (EXPERIMENTS == 1)
    {
        std::string results_path_folder = "./../results/opt_differences/";
        std::string results_path = results_path_folder + DATASET + ".txt";
        calculate_opts_differences(dataset_path, results_path, filenames);
    }
}