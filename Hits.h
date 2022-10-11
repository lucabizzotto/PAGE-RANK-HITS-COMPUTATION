#include <iostream>
#include <fstream>
#include <string>
#include <regex>
#include <map>
#include <sys/mman.h>
#include <cstdlib>
#include <cmath>
#include <algorithm>

typedef std::pair<unsigned int, unsigned int> pair;
//               val    , col position
typedef std::pair<double, unsigned int> pair_t;
//             node row value, cardinality of out links
typedef std::map<unsigned int, unsigned int> map;


class Hits
{


public:

    Hits(std::string name_file) : NAME_FILE(name_file), n_convergence(0)
    {
        // determine the number of nodes and edges in the graph
        this->number_of_nodes_and_edges();
        // read file and save values
        this->initialize();
        // compute L and L transpose and free memory of ptr pointer
        this->compute_L_and_L_t();
        // initialize hub and autority at time zero
        this->initialize_ak_hk();
    }

    // number of interations to converge
    unsigned int n_convergence;
    //autoritative score time k
    std::vector<double> a_k;
    //hub score time k
    std::vector<double> h_k;

    // custom class to execute a custom sort
    class Comparator_first_element {
    public:
        bool operator()(const pair& a, const pair& b) {
            return a.first < b.first;
        }
    };


    // custom class to execute a custom sort
    class Comparator_second_element {
    public:
        bool operator()(const pair& a, const pair& b) {
            return a.second < b.second;
        }
    };


    // compute autoritative and hub
    void compute()
    {
        this->n_convergence = 0;
        std::vector<double> tmp_a_k(this->N_nodes, 0.);
        std::vector<double> tmp_h_k(this->N_nodes, 0.);
        bool first_cicle = true;
        do
        {
            // number to converge
            this->n_convergence++;
            // COMPUTE H at time 1
            // holds the position of the actual row
            unsigned int tmp_pos_row = 0;
            // holds the position where start the next row
            unsigned int next_starting_row = this->row_ptr_L[1];
            for (unsigned int i = 0; i < this->N_edges; i++)
            {
                // when the row changes
                if (i == next_starting_row)
                {
                    // change row
                    tmp_pos_row++;
                    //                                we don't go in segmentation because the row_ptr has an extra element at the end
                    next_starting_row = this->row_ptr_L[tmp_pos_row + 1];
                }
                if (first_cicle)
                {
                    // compute value
                    tmp_h_k[this->row_ptr_not_empty_L[tmp_pos_row]] += 1;
                }
                else
                {
                    // compute value
                    tmp_h_k[this->row_ptr_not_empty_L[tmp_pos_row]] += this->a_k[this->L_ptr[i]];
                }
            }

            // COMPUTE A at time 1
            // holds the position of the actual row
            tmp_pos_row = 0;
            // holds the position where start the next row
            next_starting_row = this->row_ptr_L_t[1];
            for (unsigned int i = 0; i < this->N_edges; i++)
            {
                // when the row changes
                if (i == next_starting_row)
                {
                    // change row
                    tmp_pos_row++;
                    //                                we don't go in segmentation because the row_ptr has an extra element at the end
                    next_starting_row = this->row_ptr_L_t[tmp_pos_row + 1];
                }
                if (first_cicle)
                {
                    // compute value
                    tmp_a_k[this->row_ptr_not_empty_L_t[tmp_pos_row]] += 1;
                }
                else
                {
                    // compute value
                    tmp_a_k[this->row_ptr_not_empty_L_t[tmp_pos_row]] += this->h_k[this->L_t_ptr[i]];
                }
            }
            if (first_cicle)
            {
                first_cicle = false;
            }
            this->normalize_ak_hk(tmp_a_k, tmp_h_k);
        } while (this->distance(tmp_a_k, tmp_h_k));
    }


private:

    // holds the path of file
    const std::string NAME_FILE;
    // pointer to permanent memory
    pair* ptr;
    // holds L storing the column position
    unsigned int* L_ptr;
    // holds the starting of new line
    std::vector<unsigned int> row_ptr_L;
    // holds the empty line in L
    std::vector<unsigned int> row_ptr_not_empty_L;
    // holds L transpose storing the column position
    unsigned int* L_t_ptr;
    // holds the starting of new line
    std::vector<unsigned int> row_ptr_L_t;
    // holds the empty line in L transpose
    std::vector<unsigned int> row_ptr_not_empty_L_t;
    // holds the dimension of ptr and the number of edges
    unsigned int N_edges;
    // holds the number of nodes
    unsigned int N_nodes;
    // holds the value of min node
    unsigned int Min_node;
    // holds the value of max node
    unsigned int Max_node;


    // check the case if there is an empty row
    // we have to do because otherwise the computation of autoritative and hub aren't right in case we have at least one row complitely empty
    unsigned int check_empty_row(std::vector<unsigned int> & empty_vector, unsigned int tmp_pos)
    {
        std::vector<unsigned int>::iterator it = std::find(empty_vector.begin(), empty_vector.end(), tmp_pos);
        while (it != empty_vector.end())
        {
            if (tmp_pos == *it)
            {
                tmp_pos++;
            }
            else
            {
                break;
            }
            it++;
        }
        return tmp_pos;
    }


    // return the number of both nodes and edges in the graph
    void number_of_nodes_and_edges()
    {
        std::ifstream file = this->open_file();
        // create string to read file
        std::string line;
        // create string to store all info of file
        std::string info_string;
        // read file
        while (std::getline(file, line))
        {   // in this way we read only the description of the file
            if (line[0] != '#')
            {
                break;
            }
            else
            {
                // append all lines of description in a single string
                info_string.append(line);
            }
        }
        // close the file
        this->close_file(file);
        // create rule for regexp to find edges
        std::regex r_edges("(Edges: [0-9]*)");
        // create rule for regexp to find nodes
        std::regex r_nodes("(Nodes: [0-9]*)");
        std::regex r_min_node("(Min_node: [0-9]*)");
        std::regex r_max_node("(Max_node: [0-9]*)");
        // find the numbers
        this->N_edges = regex_exec(r_edges, info_string);
        this->N_nodes = regex_exec(r_nodes, info_string);
        this->Min_node = regex_exec(r_min_node, info_string);
        this->Max_node = regex_exec(r_max_node, info_string);
    }


    // open file
    inline std::ifstream open_file()
    {
        // create instance of file
        std::ifstream file;
        // open file
        file.open(this->NAME_FILE);
        // control the validity of file
        if (!file.is_open())
        {
            throw std::runtime_error("Could not open file");
        }
        return file;
    }


    // close file
    inline void close_file(std::ifstream& file)
    {
        file.close();
    }


    // return the value associate with regex and string
    inline int regex_exec(std::regex r, std::string info_string)
    {
        // array for regexp result
        std::smatch m;
        // compute the regexp on text
        std::regex_search(info_string, m, r);
        // retrieve the result of regexp
        std::string result = m.str(0);
        // split the result in strings with a specific break in this way we can take only the number and transform it in an integer
        std::istringstream result_line(result);
        std::vector<std::string> v_result;
        for (std::string piece; std::getline(result_line, piece, ' ');)
        {
            v_result.push_back(piece);
        }
        // return the number
        return std::stoi(v_result[1]);
    }


    // read the file and save in permanent memory
    void initialize()
    {
        // allocate the right amount of space in permanent memory to save the edges of graph as a set of pairs
        this->ptr = (pair*)mmap(NULL, this->N_edges * sizeof(pair), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, 0, 0);
        // control if the allocate was done in right way
        if (ptr == MAP_FAILED)
        {
            throw std::runtime_error("Mapping Failed\n");
        }
        // create instance of file
        std::ifstream file = this->open_file();
        // create string to read each row of file
        std::string line;
        // take the position of the ptr to insert the pair
        unsigned int i = 0;
        // read file
        while (std::getline(file, line))
        {   // in this way we don't read the description of the file
            if (line[0] != '#')
            {
                // store the pair that define the edge
                this->ptr[i] = this->take_edge_from_line(line);
                // increment the pointer
                i++;
            }
        }
        // close file
        this->close_file(file);
    }


    // from the line return the pair that compose the edge
    inline pair take_edge_from_line(std::string& line)
    {
        // split the result in strings with a specific break
        std::istringstream result_line(line);
        // contains the values for the pair
        std::vector<int> v_result;
        for (std::string piece; std::getline(result_line, piece, '\t');)
        {
            v_result.push_back(std::stoi(piece));
        }
        return pair(v_result[0], v_result[1]);
    }


    // create L
    void create_L()
    {
        // allocate the right amount of space in permanent memory to save the edges of graph as a set of pairs
        this->L_ptr = (unsigned int*)mmap(NULL, this->N_edges * sizeof(unsigned int), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, 0, 0);
        // control if the allocate was done in right way
        if (L_ptr == MAP_FAILED)
        {
            throw std::runtime_error("Mapping Failed\n");
        }
        // starting of first row
        this->row_ptr_L.push_back(0);
        this->row_ptr_not_empty_L.push_back(this->ptr[0].first - this->Min_node);
        // compute L
        for (unsigned int i = 0; i < this->N_edges; i++)
        {
            // in case of changing row
            if (i > 0 && this->ptr[i - 1].first != this->ptr[i].first)
            {
                // save the position
                row_ptr_L.push_back(i);
                this->row_ptr_not_empty_L.push_back(this->ptr[i].first - this->Min_node);
            }
            this->L_ptr[i] = this->ptr[i].second;
        }
        this->row_ptr_L.push_back(this->N_edges);
    }


    // create L transpose
    void create_L_t()
    {
        // allocate the right amount of space in permanent memory to save the edges of graph as a set of pairs
        this->L_t_ptr = (unsigned int*)mmap(NULL, this->N_edges * sizeof(unsigned int), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, 0, 0);
        // control if the allocate was done in right way
        if (L_t_ptr == MAP_FAILED)
        {
            throw std::runtime_error("Mapping Failed\n");
        }
        // starting of first row
        this->row_ptr_L_t.push_back(0);
        this->row_ptr_not_empty_L_t.push_back(this->ptr[0].second - this->Min_node);
        // compute L transpose
        for (unsigned int i = 0; i < this->N_edges; i++)
        {
            // in case of changing row
            if (i > 0 && this->ptr[i - 1].second != this->ptr[i].second)
            {
                // save the position
                row_ptr_L_t.push_back(i);
                this->row_ptr_not_empty_L_t.push_back(this->ptr[i].second - this->Min_node);
            }
            this->L_t_ptr[i] = this->ptr[i].first;
        }
        this->row_ptr_L_t.push_back(this->N_edges);
    }


    // compute matrix L and matrix L transpose
    void compute_L_and_L_t()
    {
        //sort ptr with a stable sort for the first column, in this way we are able to compute the affinity matrix
        std::stable_sort(this->ptr, this->ptr + this->N_edges, Comparator_first_element());
        this->create_L();
        //sort again by second element 
        std::stable_sort(this->ptr, this->ptr + this->N_edges, Comparator_second_element());
        this->create_L_t();
        // free memory 
        if (munmap(this->ptr, this->N_edges) != 0)
        {
            throw std::runtime_error("Free memory failed\n");
        }
    }

    
    //initalize hub score and autoritative score
    void initialize_ak_hk()
    {
        this->a_k.clear();
        this->a_k.resize(this->N_nodes, 0.);
        this->h_k.clear();
        this->h_k.resize(this->N_nodes, 0.);
    }


    //to normalize the vector to have a probability distribution
    void normalize_ak_hk(std::vector<double>& vec_to_normalize_ak, std::vector<double>& vec_to_normalize_hk)
    {
        double sum_a_k = 0.0;
        double sum_h_k = 0.0;
        for (int i = 0; i < this->N_nodes; i++)
        {
            sum_a_k += vec_to_normalize_ak[i];
            sum_h_k += vec_to_normalize_hk[i];
        }
        for (int i = 0; i < this->N_nodes; i++)
        {
            vec_to_normalize_ak[i] = vec_to_normalize_ak[i] / sum_a_k;
            vec_to_normalize_hk[i] = vec_to_normalize_hk[i] / sum_h_k;
        }
    }


    // compute the distance between autoritative and hub at time t and t+1
    bool distance(std::vector<double>& ak_tmp, std::vector<double>& hk_tmp)
    {
        double distance_ak = 0.;
        double distance_hk = 0.;
        double mean = 0.0;

        for (int i = 0; i < this->N_nodes; i++)
        {
            distance_ak = distance_ak +  std::abs(this->a_k[i] - ak_tmp[i]);
            distance_hk = distance_hk +  std::abs( this->h_k[i] - hk_tmp[i]);
        }
        mean = (distance_ak + distance_hk) / 2.;
        //I have to keep iterating
        if (mean > std::pow(10, -10))
        {
            this->a_k = ak_tmp;
            this->h_k = hk_tmp;
            //clean the tmp for next iteration
            ak_tmp.clear();
            ak_tmp.resize(this->N_nodes, 0.);
            hk_tmp.clear();
            hk_tmp.resize(this->N_nodes, 0.);
            return true;
        }
        return false;
    }


};