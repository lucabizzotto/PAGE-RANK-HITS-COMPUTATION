#include <iostream>
#include <fstream>
#include <string>
#include <regex>
#include <map>
#include <sys/mman.h>
#include <cstdlib>
#include <cmath>

typedef std::pair<unsigned int, unsigned int> pair;
//               val    , col position
typedef std::pair<double, unsigned int> pair_t;
//             node row value, cardinality of out links
typedef std::map<unsigned int, unsigned int> map;


class PageRank
{
    public:

        PageRank(std::string name_file, double d) : NAME_FILE(name_file), TELEPORTATION_PROB(d), n_convergence(0)
        {
            // determine the number of nodes and edges in the graph
            number_of_nodes_and_edges();
            // read file and save values
            initialize();
            // create cardinality map and dangling vector
            create_dangling_and_cardinality_map();
            // build the transpose matrix
            create_transpose_matrix();
            // this is the prestige at time zero
            initialize_p_k();
        }


        // vector of prestige
        std::vector<double> p_k;
        // number of interations to converge
        unsigned int n_convergence;


        // custom class to execute a custom sort
        class Compararator_first_element {
        public:
            bool operator()(const pair& a, const pair& b) {
                return a.first < b.first;
            }
        };


        // custom class to execute a custom sort
        class Compararator_second_element {
            public:
                bool operator()(const pair& a, const pair& b) {
                    return a.second < b.second;
                }
        };


        // compute pagerank
        void compute()
        {
            this->n_convergence = 0;
            // vector that indicates the moltiplication between A transpose and p_k initialized with zeros
            std::vector<double> tmp_p_k(this->N_nodes, 0.);
            do
            {
                // number to converge
                this->n_convergence++;
                // COMPUTE A^t * P_k
                // computing the moltiplication between A transpose and p_k
                // holds the position of the actual row
                unsigned int tmp_pos_row = 0;
                // holds the position where start the next row
                unsigned int next_starting_row = this->row_ptr[1];
                for (unsigned int i = 0; i < this->N_edges; i++)
                {
                    // when the row changes
                    if (i == next_starting_row)
                    {
                        // change row
                        tmp_pos_row++;
                        //                                we don't go in segmentation because the row_ptr has an extra element at the end
                        next_starting_row = this->row_ptr[tmp_pos_row + 1];
                    }
                    // compute value
                    tmp_p_k[this->row_ptr_transpose_not_empty[tmp_pos_row]] += this->ptr_transpose[i].first * this->p_k[this->ptr_transpose[i].second - this->Min_node];
                }
                // COMPUTE DANGLING * P_k * 1/N
                // indicates the moltiplication between dangling and p_k initialized with zeros (we need only one number because each row has the same value
                double D_x_pk = 0.;
                for (auto it : this->dangling_nodes)
                {
                    D_x_pk += this->p_k[it - this->Min_node] * (1. / this->N_nodes);
                }
                // COMPUTE P_k at time t + 1
                for (std::size_t i = 0; i < tmp_p_k.size(); i++)
                {
                    //            (         (A^t + Dangling) * P_k) * d           ) + (( 1    -    d) e/n                           )   
                    tmp_p_k[i] = ((D_x_pk + tmp_p_k[i]) * this->TELEPORTATION_PROB) + (1 - this->TELEPORTATION_PROB) / this->N_nodes;
                }
            } while (distance_between_p(tmp_p_k));
        }

        
    private:

        // holds the path of file
        const std::string NAME_FILE;
        // value of probability to not use teleportation
        const double TELEPORTATION_PROB;
        // pointer to permanent memory
        pair *ptr;
        // pointer to permanent memory
        pair_t* ptr_transpose;
        // holds the dimension of ptr and the number of edges
        unsigned int N_edges;
        // holds the number of nodes
        unsigned int N_nodes;
        // holds the value of min node
        unsigned int Min_node;
        // holds the value of max node
        unsigned int Max_node;
        // map to maintain cardinality of out links
        map cardinality_map;
        // vector that contains the first position in transpose matrix about the element 
        std::vector<unsigned int> row_ptr;
        // vector of dangling nodes
        std::vector<unsigned int> dangling_nodes;
        // holds the empty line in L transpose
        std::vector<unsigned int> row_ptr_transpose_not_empty;

        // initialize p_k at time 0
        inline void initialize_p_k()
        {
            this->p_k.clear();
            this->p_k.resize(this->N_nodes, 1. / this->N_nodes);
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


        // return the number of both nodes and edges in the graph
        void number_of_nodes_and_edges()
        {   
            std::ifstream file = this->open_file();
            // create string to read file
            std::string line;
            // create string to store all info of file
            std::string info_string;
            // read file
            while(std::getline(file, line))
            {   // in this way we read only the description of the file
                if(line[0] != '#')
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


        // insert dangling in right way
        inline void insert_dangling(unsigned int start, unsigned int end)
        {
            // we insert the nodes that are in the middle
            for (unsigned int u = start + 1; u <= end; u++)
            {
                this->dangling_nodes.push_back(u);
            }
        }


        // read the file and save in permanent memory
        void initialize()
        {
            // allocate the right amount of space in permanent memory to save the edges of graph as a set of pairs
            this->ptr = (pair*)mmap (NULL, this->N_edges*sizeof(pair), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, 0, 0);
            // control if the allocate was done in right way
            if(ptr == MAP_FAILED)
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
            while(std::getline(file, line))
            {   // in this way we don't read the description of the file
                if(line[0] != '#')
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


        // create the dangling vector and cardinality map
        void create_dangling_and_cardinality_map()
        {
            // sort ptr with a stable sort for the first column, in this way we are able to find dangling and cardinality for each node
            std::stable_sort(this->ptr, this->ptr + this->N_edges, Compararator_first_element());
            // maintain the value of node to store the amount of out link that each node has
            int node = -1;
            // maintain the cardinality of node, so the degree of node
            unsigned int cardinality_of_node = 0;
            // counter to know when there is a dangling nodes
            int predecessor = -1;
            // read all ptr ordered
            for (unsigned int i = 0; i < this->N_edges; i++)
            {
                // the first case that the node isn't initialized
                if (node == -1)
                {
                    // set to 1 because we have just seen an edge of the current node
                    cardinality_of_node = 1;
                    node = this->ptr[i].first;
                    predecessor = this->ptr[i].first;
                }
                // means that we are still saying the same node so we only increase its degree
                else if (this->ptr[i].first == static_cast<unsigned int>(node))
                {
                    cardinality_of_node++;
                }
                // means that the node is changed
                else
                {
                    this->cardinality_map[static_cast<unsigned int> (node)] = cardinality_of_node;
                    // set to 1 because we have just seen an edge of the current node
                    cardinality_of_node = 1;
                    node = this->ptr[i].first;
                    predecessor = this->ptr[i - 1].first;
                }
                // check the presence of dangling nodes
                if (this->ptr[i].first - predecessor > 1 && predecessor != -1)
                {
                    insert_dangling(predecessor, this->ptr[i].first - 1);
                    // to avoid double check of the same node
                    predecessor = this->ptr[i].first;
                }
            }
            // to add the last danglings  
            if ((this->Max_node) - (this->ptr[this->N_edges - 1].first) > 0)
            {
                insert_dangling(this->ptr[this->N_edges - 1].first, this->Max_node);
            }
            // to add the last element in cardinality map
            this->cardinality_map[node] = cardinality_of_node;
        }


        // build the transpose matrix
        void create_transpose_matrix()
        {
            // sort ptr with a stable sort for the second column, in this way we are able to compute the transpose matrix
            std::stable_sort(this->ptr, this->ptr + this->N_edges, Compararator_second_element());
            // allocate the right amount of space in permanent memory to save the transpose matrix 
            this->ptr_transpose = (pair_t*)mmap(NULL, this->N_edges * sizeof(pair_t), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, 0, 0);
            // starting of first row
            this->row_ptr.push_back(0);
            this->row_ptr_transpose_not_empty.push_back(this->ptr[0].second - this->Min_node);
            // populate the transpose matrix
            for (unsigned int i = 0; i < this->N_edges; i++) {
                // in case of changing row
                if (i > 0 && this->ptr[i - 1].second != this->ptr[i].second)
                {
                    // save the position
                    row_ptr.push_back(i);
                    this->row_ptr_transpose_not_empty.push_back(this->ptr[i].second - this->Min_node);
                }
                // store value 1/k
                this->ptr_transpose[i] = pair_t(1. / this->cardinality_map[this->ptr[i].first], this->ptr[i].first);
            }
            // store the value that indicates the end of matrix
            this->row_ptr.push_back(this->N_edges);
            // free memory 
            if (munmap(this->ptr, this->N_edges) != 0)
            {
                throw std::runtime_error("Free memory failed\n");
            }
        }


        // calculate the distance between p at time t and time t+1
        inline bool distance_between_p( std::vector<double> & tmp_p_k)
        {
            double distance = 0.;
            for (unsigned int i = 0; i < tmp_p_k.size(); i++)
            {
                distance += std::abs(tmp_p_k[i] - this->p_k[i]); // norm 1
            }
            // assign to p_k the p_k+1
            this->p_k = tmp_p_k;
            // initialize vector
            tmp_p_k.clear();
            tmp_p_k.resize(this->N_nodes, 0.);
            return distance > std::pow(10, -10);
        }
};