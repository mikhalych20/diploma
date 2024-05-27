
#include <vector>
#include <iostream>
#include <fstream>
#include <queue>

//fstream file("output.txt");
//fstream file("output.txt");
using namespace std;
ifstream file("C:\\Users\\mimik\\source\\repos\\generate_data\\generate_data\\output.txt");
ofstream out_file("C:\\Users\\mimik\\source\\repos\\generate_data\\generate_data\\intput_1.txt");
class graph
{
public:
	int size;
	vector<vector<int>> matrix_graph;
	vector<int> points;
	graph()
	{

	}
	graph(int size_, int point)
	{
		size = size_;

	}
	graph(int size_)
	{
		size = size_;
		vector<int> tmp_points(size);
		for (int i = 0; i < size; ++i)
		{

			points.push_back(i);
			vector<int> vect(size);
			matrix_graph.push_back(vect);
			for (int j = 0; j < size; ++j)
			{
				int tmp;
				file >> tmp;
				matrix_graph[i][j] =tmp;
			}
		}
	}

	void print_graph()

	{
		out_file << "POINTS ";
		for (int i = 0; i < this->size; ++i)
		{
			out_file << this->points[i];
		}
		out_file << "\nGRAPH\n";
		for (int i = 0; i < this->size; ++i)
		{
			for (int j = 0; j < this->size; ++j)
			{
				out_file << this->matrix_graph[i][j];

			}
			out_file << "\n";
		}
		out_file << "\n";
	}


};

bool check_C4(vector<vector<int>>& gr)
{
	bool res = true;
	for (int i = 0; i < 4; i++)
	{
		int tmp = 0;
		for (int j = 0; j < 4; j++)
		{
			if (gr[i][j] == 1)
				tmp++;
		}
		if (tmp != 2)
		{
			res = false;
			break;
		}
	}
	return res;

}

bool check_P4(vector<vector<int>>& gr)
{
	bool res = true;
	int tmp_1 = 0;
	int tmp_2 = 0;
	for (int i = 0; i < 4; i++)
	{
		int tmp = 0;

		for (int j = 0; j < 4; j++)
		{
			if (gr[i][j] == 1)
				tmp++;
		}
		if (tmp == 1)
		{
			tmp_1++;
		}
		if (tmp == 2)
		{
			tmp_2++;
		}
	}
	if ((tmp_1 != 2) || (tmp_2 != 2))
		res = false;
	return res;
}

bool check_2K2(vector<vector<int>>& gr)
{
	bool res = true;
	for (int i = 0; i < 4; i++)
	{
		int tmp = 0;
		for (int j = 0; j < 4; j++)
		{
			if (gr[i][j] == 1)
				tmp++;
		}
		if (tmp != 1)
		{
			res = false;
			break;
		}
	}
	return res;
}

bool check_for_threshold_graph_with_subgraphs(graph& checked_graph)
{
	bool result = true;
	vector<int> tmp(4, 0);
	vector<vector<int>> subgraph(4, tmp);
	for (int i = 0; i < checked_graph.size; ++i)
	{
		for (int j = i + 1; j < checked_graph.size; ++j)
		{
			for (int k = j + 1; k < checked_graph.size; ++k)
			{
				for (int l = k + 1; l < checked_graph.size; ++l)
				{
					subgraph[0][1] = subgraph[1][0] = checked_graph.matrix_graph[i][j];
					subgraph[0][2] = subgraph[2][0] = checked_graph.matrix_graph[i][k];
					subgraph[0][3] = subgraph[3][0] = checked_graph.matrix_graph[i][l];
					subgraph[1][2] = subgraph[2][1] = checked_graph.matrix_graph[j][k];
					subgraph[1][3] = subgraph[3][1] = checked_graph.matrix_graph[j][l];
					subgraph[2][3] = subgraph[3][2] = checked_graph.matrix_graph[k][l];

					if (check_C4(subgraph) || check_2K2(subgraph) || check_P4(subgraph))
					{
						result = false;
						i = j = k = l = checked_graph.size;;
					}


				}
			}
		}
	}


	return result;
}

///////////////////////////////ESTIMATE
bool check_for_empty(graph& checked_graph)
{
	bool result = true;
	for (int i = 0; i < checked_graph.matrix_graph.size(); ++i)
		for (int j = i; j < checked_graph.matrix_graph.size(); ++j) {
			if (checked_graph.matrix_graph[i][j] == 1)
			{
				result = false;
				break;
			}
		}

	return result;
}

struct est_node
{
	graph data;
	est_node* root;
	est_node* left_son;
	est_node* right_son;
	

	
};

est_node get_est_root(graph& getted_graph)
{
	est_node result;
	result.data = getted_graph;
	result.root = result.left_son = result.right_son = NULL;
	return result;
}
void add_est_node(est_node* node, est_node* root, int type)
{
	if (type == 0)
	{
		root->left_son = node;
		node->root = root;
	}
	if (type == 1)
	{
		root->right_son = node;
		node->root = root;
	}
}
graph right_decompozition(graph& current_graph, int pozition)
{
	graph result(current_graph.size, 1);
	//vector<int> tmp(current_graph.matrix_graph.size() - 1, 0);
	vector<vector<int>> first_m = current_graph.matrix_graph;
	for (int j = 0; j < first_m.size(); ++j)
	{
		//	cout << "+++ " << j << "\n";
		if (first_m[j][pozition] == 1)
		{
			//cout << "--- "<<j<<"\n";
			for (int i = 0; i < first_m.size(); ++i)
			{
				first_m[i].erase(first_m[i].begin() + j);

			}
			first_m.erase(first_m.begin() + j);
			if (pozition > j)
				pozition = pozition - 1;
			j = j - 1;
		}
	}
	result.matrix_graph = first_m;
	result.size = first_m.size();
	return result;
}

graph left_decompozition(graph& current_graph, int pozition)
{
	graph first(current_graph.size, 1);
	vector<vector<int>> first_m = current_graph.matrix_graph;
	for (int i = 0; i < current_graph.matrix_graph.size(); ++i)
	{
		first_m[i].erase(first_m[i].begin() + pozition);

	}
	first_m.erase(first_m.begin() + pozition);

	first.matrix_graph = first_m;
	return first;
}

int est_tree(est_node* current_node, int current_est_number)
{
	int result_number = current_est_number;
	if (check_for_empty(current_node->data))
	{
		result_number = current_node->data.matrix_graph.size();
	}
	else
	{
		/////���� ��������������� �������
		int find_i = 0;
		for (int i = 0; i < current_node->data.matrix_graph.size(); ++i)
		{
			for (int j = i; j < current_node->data.matrix_graph.size(); ++j)
				if (current_node->data.matrix_graph[i][j] == 1)
				{
					find_i = i;
					i = j = current_node->data.matrix_graph.size();
				}
		}
		//////////////
		//// ������ ������ ����
		graph left_son = left_decompozition(current_node->data, find_i);
		est_node left_son_ = get_est_root(left_son);
		add_est_node(&left_son_, current_node, 0);
		/////////
		/////��������� �� �� ����� � ����� ����� � ��������
		int left_result_number = est_tree(&left_son_, result_number);

		/////////////////
		/////������ ������� ����
		graph right_son = right_decompozition(current_node->data, find_i);
		est_node right_son_ = get_est_root(right_son);
		add_est_node(&right_son_, current_node, 1);
		///////////////
		///// ��������� �� �� ����� ��� ������� ����
		int right_result_number = est_tree(&right_son_, result_number);
		result_number = max(left_result_number, right_result_number);
		////////////
	}
	if (result_number < current_est_number)
	{
		result_number = current_est_number;
	}
	return result_number;
}
int get_estimate_threshold_number(graph& checked_graph)
{
	int result = 0;
	if (check_for_empty(checked_graph))
		result = checked_graph.size;
	else
	{
		est_node root = get_est_root(checked_graph);
		result = est_tree(&root, 1);
	}
	/////////////////////////////////�������� ���������� ������!!!
	result = checked_graph.size - result;
	//////////////////////
	return result;
}


///////////////////////////////////////////////////////////REAL_NUMBER

struct Node
{
	vector<graph> data;
	vector<Node*> leafs;
	Node* root;
	int count_of_trashold;
	int layer;
};

Node create_node(vector<graph>& data)
{
	Node root_;
	root_.root = NULL;
	root_.data = data;
	vector<Node*> leafs;
	root_.leafs = leafs;
	root_.count_of_trashold = 0;
	return root_;
}

void add_Node(Node* node, Node* root)
{
	node->root = root;
	root->leafs.push_back(node);
	node->layer = root->layer + 1;
}

Node get_root(graph& getted_graph)
{
	Node root_;
	root_.layer = 0;
	root_.root = NULL;
	vector<graph> data;
	data.push_back(getted_graph);
	root_.data = data;
	vector<Node*> leafs;
	root_.leafs = leafs;
	root_.count_of_trashold = 0;
	return root_;


}

vector<graph>  decompozition(graph& current_graph, int first_deleted_index, int second_deleted_index)
{
	/*cout << "DECOMPOZITIO\n";

	print_graph(current_graph);

	cout << "DECOMPOZITIO\n";*/
	vector<graph> result(2);
	vector<int> first_points = current_graph.points;
	vector<int> second_points = current_graph.points;
	graph first(current_graph.matrix_graph.size() - 1, 1);
	graph second(current_graph.matrix_graph.size() - 1, 1);
	//vector<vector<int>> first_m = current_graph.matrix_graph;
	first.matrix_graph = current_graph.matrix_graph;
	//vector<vector<int>> second_m = current_graph.matrix_graph;
	second.matrix_graph = current_graph.matrix_graph;
	for (int i = 0; i < current_graph.matrix_graph.size(); ++i)
	{
		first.matrix_graph[i].erase(first.matrix_graph[i].begin() + first_deleted_index);

	}
	first.matrix_graph.erase(first.matrix_graph.begin() + first_deleted_index);
	first_points.erase(first_points.begin() + first_deleted_index);

	for (int i = 0; i < current_graph.matrix_graph.size(); ++i)
	{
		second.matrix_graph[i].erase(second.matrix_graph[i].begin() + second_deleted_index);

	}
	second.matrix_graph.erase(second.matrix_graph.begin() + second_deleted_index);
	second_points.erase(second_points.begin() + second_deleted_index);


	/*for (int i = 0; i < current_graph.matrix_graph.size() - 1; ++i)
		for (int j = 0; j < current_graph.matrix_graph.size() - 1; ++j)
			cin >> first_m[i][j];*/

			//first_m[0][1] = first_m[0][2] = first_m[0][3] = first_m[3][0] = first_m[2][0] = first_m[1][0] = 1;
			//vector<vector<int>> second_m(current_graph.matrix_graph.size() - 1, tmp);
			/*for (int i = 0; i < current_graph.matrix_graph.size() - 1; ++i)
				for (int j = 0; j < current_graph.matrix_graph.size() - 1; ++j)
					cin >> second_m[i][j];*/
	//first.matrix_graph = first_m;
	//second_m[0][1] = second_m[1][0] = second_m[2][1] = second_m[1][2] = 1;
	//second.matrix_graph = second_m;

	first.points = first_points;
	second.points = second_points;


	////���� �������
	////����� ���� �� 2
	result[0] = first;
	result[1] = second;
	return result;
}
/*vector<graph> delete_dublicate(vector<graph> data)
{
	vector<graph> new_data = data;
	for (int i = 0; i < data.size(); ++i)
	{
		for (int j = i + 1; j < new_data.size(); ++j)
		{
			if (new_data[j] == data[i])
			{
				new_data.erase(new_data.begin() + j);
				j = j - 1;
			}
		}
	}
	return new_data;

}*/

bool find(vector<int>& d, int value)
{
	/*vector<int> data = d;
	int val = value;
	bool result = true;
	int l = 0;
	int r = data.size();
	int res;
	while (r > l) {
		int m = (l + r) / 2;    //������������� �������!

		if (data[m] < val) {
			l = m + 1;
		}
		else if (data[m] > val) {
			r = m - 1;
		}
		else {
			//res = m;
			return true;
		}
	}
	//cout << data.size() << " " << l<<"\n";
	if (l < data.size())
	{
		if (data[l] == val) {
			res = l;
		}
		else {
			res = -1;
		}
		if (res == -1)
			result = false;
	}

	
*/
	bool result = false;
	for (int i = 0; i < d.size(); ++i)
	{
		if (d[i] == value)
			result = true;
	}
	return result;

}

bool check_in(vector<int>& first, vector<int>& second)
{
	//cout<<"gfdsa";
	bool result = false;
	if (first.size() < second.size())
	{
		int check = 0;
		for (int i = 0; i < first.size(); ++i)
		{
			if (find(second, first[i]))
			{
				check++;
			}
		}
		if (check == first.size())
			result = true;
	}
	else
	{
		int check = 0;
		for (int i = 0; i < second.size(); ++i)
		{
			if (find(first, second[i]))
			{
				check++;
			}
		}
		if (check == second.size())
			result = true;
	}
	return result;
}

vector<graph> delete_subgraphs (Node* current_node)
{

	//cout<<"rgfregf";
	//out_file << "SUBGRAPh\n";
	vector<graph> new_data = current_node->data;

	for (int i = 0; i < new_data.size(); ++i)
	{
	
		for (int j = i+1; j < new_data.size(); ++j)
		{
			//	cout << "bin"<< new_data.size();
			/*out_file << "\n";
			for (int k = 0; k < new_data[i].points.size(); ++k)
			{
				out_file << new_data[i].points[k];
			}
			out_file << "\n";

			for (int k = 0; k < new_data[j].points.size(); ++k)
			{
				out_file << new_data[j].points[k];
			}

			out_file << "\n";*/

			if ( (check_in(new_data[i].points, new_data[j].points)))
			{
				//out_file << "\nDELETED\n";
				if (new_data[i].points.size() > new_data[j].points.size())
				{
					if (check_for_threshold_graph_with_subgraphs(new_data[j]))
					{
						current_node->count_of_trashold--;
					}

					new_data.erase(new_data.begin() + j);

					j--;
				}
				else
				{
					if (check_for_threshold_graph_with_subgraphs(new_data[i]))
					{
						current_node->count_of_trashold--;
					}

					new_data.erase(new_data.begin() + i);

					j--;
				}
				//	cout << "bin"<< new_data.size();
			if(current_node->count_of_trashold <0)
			{
				cout<<"!!!!!!";
			}

			}
			/*else
			{
				out_file << "\nNOT DELETED\n";
			}*/
			
		}
	}
//	cout << "delete trashold" << new_data.size() << " " << current_node->count_of_trashold << "\n";

	return new_data;

}

vector<int>get_points(Node* current_node, int i)
{
	vector<int> point(4,-1);
	vector<int> tmp(4, 0);
	vector<vector<int>> subgraph(4, tmp);
	for (int ii = 0; ii < current_node->data[i].size; ++ii)
	{
		for (int j = ii + 1; j < current_node->data[i].size; ++j)
		{
			for (int k = j + 1; k < current_node->data[i].size; ++k)
			{
				for (int l = k + 1; l < current_node->data[i].size; ++l)
				{
				subgraph[0][1] = subgraph[1][0] = current_node->data[i].matrix_graph[ii][j];
					subgraph[0][2] = subgraph[2][0] = current_node->data[i].matrix_graph[ii][k];
					subgraph[0][3] = subgraph[3][0] = current_node->data[i].matrix_graph[ii][l];
					subgraph[1][2] = subgraph[2][1] = current_node->data[i].matrix_graph[j][k];
					subgraph[1][3] = subgraph[3][1] = current_node->data[i].matrix_graph[j][l];
					subgraph[2][3] = subgraph[3][2] = current_node->data[i].matrix_graph[k][l];

					if (check_C4(subgraph) || check_2K2(subgraph) || check_P4(subgraph))
					{
						point[0] = ii;
						point[1] = j;
						point[2] = k;
						point[3] = l;
						ii = j = k = l = current_node->data[i].size;
						subgraph.clear();
					}


				}
			}
		}


	}
	return point;
}









int result_tree_long(Node* current_node, int start_index_in_data, int current_number, queue<Node>& qu)
{
	cout << current_number<<" "<<qu.size()<<" "<< current_node->data.size()<<" "<< current_node->layer<<"\n";
	//qu.pop();
	//cout << current_number << qu.size() << current_node->data.size()<< "jdjzskjdn\n";

	///vector<graph> result_decompozition;
	int find_i = -1, find_j = -1, find_k = -1, find_l = -1;
	/*if (current_node->root != NULL) {
		cout << "ROOT IS\n";
		for (int k = 0; k < current_node->root->data.size(); ++k)
		{
			cout << "           " << k << "\n";
			print_graph(current_node->root->data[k]);
			cout << "         \n";
		}
		cout << "\n";
	}

	cout << "DATA IS\n";


	for (int k = 0; k < current_node->data.size(); ++k)
	{
		cout << "         "<<k<<"\n";
		print_graph(current_node->data[k]);
		cout << "           \n";
	}

	cout << "NUMBER IS " << current_number<<"\n";
	cout << "CURRENT TRASHOLD NUMBER IS " << current_node->count_of_trashold << "\n";
	*/
	vector<int> point(4, -1);
	int result_number = current_number;

	for (int i = start_index_in_data; i < current_node->data.size(); ++i)
	{
		//cout <<"DATA"<< current_node->data.size() << "\n";



		point = get_points(current_node, i);
		find_i = point[0];
		find_j = point[1];
		find_k = point[2];
		find_l = point[3];

		point.clear();
		

		if (find_i == -1)
		{
			current_node->count_of_trashold++;
			//cout << "NEW COUNT IS " << current_node->count_of_trashold << "\n";
		}
		else
		{
	//		cout << point[0] << point[1] << point[2] << point[3];
			///////////////////////////////////////////////////////////////



			///������������, �������� ����� ������ data
		//	cout << "DATA" << current_node->data.size() << "\n";
			vector<graph>tmp_data;
			vector<graph>new_data(current_node->data.size() + 1);
			for (int j = 0; j < i; ++j)
				new_data[j] = current_node->data[j];
			for (int j = i + 1; j < current_node->data.size(); ++j)
				new_data[j + 1] = current_node->data[j];


			if (current_node->data[i].matrix_graph[find_i][find_j] == 0)
			{
				tmp_data = decompozition(current_node->data[i], find_i, find_j);
				new_data[i] = tmp_data[0];
				new_data[i + 1] = tmp_data[1];

				//////////////////////////////////////////////////
				Node new_node = create_node(new_data);
				//new_node.count_of_trashold= current_node->count_of_trashold;

				add_Node(&new_node, current_node);
				qu.push(new_node);
				//cout << current_number << qu.size() << "\n";
				//cout << current_node->data.size() << "\n";
				//cout << "DATA" << current_node->data.size() << "\n";


				//result_number = result_tree(&new_node, 0, result_number);
			}

			if (current_node->data[i].matrix_graph[find_i][find_k] == 0)
			{
				tmp_data = decompozition(current_node->data[i], find_i, find_k);
				new_data[i] = tmp_data[0];
				new_data[i + 1] = tmp_data[1];

				//////////////////////////////////////////////////
				Node new_node = create_node(new_data);
				//new_node.count_of_trashold= current_node->count_of_trashold;

				add_Node(&new_node, current_node);
				qu.push(new_node);
				//cout << current_number << qu.size() << "\n";
				//cout << current_node->data.size() << "\n";
				//cout << "DATA" << current_node->data.size() << "\n";


			//	result_number = result_tree(&new_node, 0, result_number);
			}

			if (current_node->data[i].matrix_graph[find_i][find_l] == 0)
			{
				tmp_data = decompozition(current_node->data[i], find_i, find_l);
				new_data[i] = tmp_data[0];
				new_data[i + 1] = tmp_data[1];

				//////////////////////////////////////////////////
				Node new_node = create_node(new_data);
				//new_node.count_of_trashold= current_node->count_of_trashold;

				add_Node(&new_node, current_node);
				qu.push(new_node);
				//cout << current_number << qu.size() << "\n";
				//cout << current_node->data.size() << "\n";
				//result_number = result_tree(&new_node, 0, result_number);
				//cout << "DATA" << current_node->data.size() << "\n";

			}

			if (current_node->data[i].matrix_graph[find_j][find_k] == 0)
			{
				tmp_data = decompozition(current_node->data[i], find_j, find_k);
				new_data[i] = tmp_data[0];
				new_data[i + 1] = tmp_data[1];

				//////////////////////////////////////////////////
				Node new_node = create_node(new_data);
				//new_node.count_of_trashold= current_node->count_of_trashold;

				add_Node(&new_node, current_node);
				qu.push(new_node);
				//cout << current_number << qu.size() << "\n";
				//cout << current_node->data.size() << "\n";
				//cout << "DATA" << current_node->data.size() << "\n";


				//result_number = result_tree(&new_node, 0, result_number);
			}

			if (current_node->data[i].matrix_graph[find_j][find_l] == 0)
			{
				tmp_data = decompozition(current_node->data[i], find_j, find_l);
				new_data[i] = tmp_data[0];
				new_data[i + 1] = tmp_data[1];

				//////////////////////////////////////////////////
				Node new_node = create_node(new_data);
				//new_node.count_of_trashold= current_node->count_of_trashold;

				add_Node(&new_node, current_node);
				qu.push(new_node);
				//cout << current_number << qu.size() << "\n";
				//cout << current_node->data.size() << "\n";
				//cout << "DATA" << current_node->data.size() << "\n";


				//result_number = result_tree(&new_node, 0, result_number);
			}

			if (current_node->data[i].matrix_graph[find_k][find_l] == 0)
			{
				tmp_data = decompozition(current_node->data[i], find_k, find_l);
				new_data[i] = tmp_data[0];
				new_data[i + 1] = tmp_data[1];

				//////////////////////////////////////////////////
				Node new_node = create_node(new_data);
				//new_node.count_of_trashold= current_node->count_of_trashold;

				add_Node(&new_node, current_node);
				qu.push(new_node);
				//cout << current_number << qu.size() << "\n";
				//cout << current_node->data.size() << "\n";
				//cout << "DATA" << current_node->data.size() << "\n";

				//result_number = result_tree(&new_node, 0, result_number);
			}

			///��������� ����� ������� � ������ � �������, count_of_trashold �������
			/// ��������� ��� �� ������� ��� ����� ������� � start_index = i
		//	cout << "DATA" << current_node->data.size() << "\n";

		}
		//cout << "DATA" << current_node->data.size() << "\n";

	}
	//vector<graph> new_data= delete_dublicate(current_node->data);
	//current_node->data = delete_subgraphs(current_node);
//	cout << "DATA" << current_node->data.size() <<" "<< current_node->count_of_trashold<< "\n";
	if (current_node->count_of_trashold == current_node->data.size())
	{
		cout << "1111\n";

		//cout << "PPPPPPPP"<< current_node->count_of_trashold<< current_node->data.size();

		//->data = delete_subgraphs(current_node);
	 /* if(current_node->count_of_trashold == 2)*/
		

		//current_node->data = delete_subgraphs(current_node);
		// cout<<"THIS";
		 //current_node->data = delete_dublicate(current_node->data);
		 //
		if (result_number >= current_node->data.size())
		{
			result_number = current_node->data.size();
			//cout << result_number;
			//int layer = current_node->layer;
			//if (!(qu.empty()))
			//   qu.pop();
			//if(!(qu.empty()))
			//	if(qu.front().layer>layer)
					while (!qu.empty()) {
						qu.pop();
					}

			//cout << "\nNUMBER IS CHANGE\n"<< result_number<<"\n";

			/*for (int ind_in_data = 0; ind_in_data < current_node->data.size(); ind_in_data++)
			{
				cout << "index " << ind_in_data << "\n";
				print_graph(current_node->data[ind_in_data]);
			}*/
			/*for (int ind_in_data = 0; ind_in_data < current_node->data.size(); ind_in_data++)
			{
				out_file << "index " << ind_in_data << "\n";
				current_node->data[ind_in_data].print_graph();
			}*/

			/*while (!qu.empty()) {
				qu.pop();
			}*/
		}

	}
	else
	{
		if (!(qu.empty()))
			qu.pop();
	}
	while (!(qu.empty()))
	{

		result_number = result_tree_long(&qu.front(), 0, result_number, qu);
	}

	return result_number;
}



int get_threshold_number_with_decition_tree_long(graph& checked_graph)
{
	int result = 0;
	bool are_trashhold = false;
	Node root = get_root(checked_graph);
	queue<Node> qu;
	if (check_for_threshold_graph_with_subgraphs(checked_graph))
	{
		result = 1;
	}
	else
	{
		//		vector<graph> result_decompozition;
		qu.push(root);
		result = result_tree_long(&root, 0, checked_graph.size, qu);
		//		result = size(result_decompozition);
	}
	return result;
}



int get_threshold_number_with_decition_tree_long_cycle(graph& checked_graph)
{
	//cout<<"thiss";
	int result= checked_graph.size-1;
	Node root = get_root(checked_graph);
	queue<Node&> qu;
	qu.push(root);
	Node current_node;
	int find =0;
	while(!(qu.empty()))
	{
			int find_i = -1, find_j = -1, find_k = -1, find_l = -1;
			current_node = qu.front();
			vector<int> point(4, -1);
			//cout << result<<" "<<qu.size()<<" "<< current_node.data.size()<<" "<< current_node.layer<<"\n";
			for(int i =0; i < current_node.data.size();++i)
			{
				point = get_points(&current_node, i);
				find_i = point[0];
				find_j = point[1];
				find_k = point[2];
				find_l = point[3];
				point.clear();
				if (find_i == -1)
				{
					current_node.count_of_trashold++;
					//cout << "NEW COUNT IS " << current_node.count_of_trashold << "\n";
				}
			}
			current_node.data = delete_subgraphs(&current_node);
			//cout << result<<" "<<qu.size()<<" "<< current_node.data.size()<<" "<< current_node.count_of_trashold<< " " <<find<<"\n";

			for(int i =0; i < current_node.data.size();++i)
			{
				point = get_points(&current_node, i);
				find_i = point[0];
				find_j = point[1];
				find_k = point[2];
				find_l = point[3];
				point.clear();
				if (find_i == -1)
				{
					//current_node.count_of_trashold++;
					//cout << "NEW COUNT IS " << current_node.count_of_trashold << "\n";
				}
				else
				{
					vector<graph>tmp_data;
					vector<graph>new_data(current_node.data.size() + 1);
					for (int j = 0; j < i; ++j)
						new_data[j] = current_node.data[j];
					for (int j = i + 1; j < current_node.data.size(); ++j)
						new_data[j + 1] = current_node.data[j];

						if (current_node.data[i].matrix_graph[find_i][find_j] == 0)
					{
						tmp_data = decompozition(current_node.data[i], find_i, find_j);
						new_data[i] = tmp_data[0];
						new_data[i + 1] = tmp_data[1];

						//////////////////////////////////////////////////
						Node new_node = create_node(new_data);
						//new_node.count_of_trashold= current_node->count_of_trashold;

						add_Node(&new_node, &current_node);
						qu.push(new_node);
						//cout << current_number << qu.size() << "\n";
						//cout << current_node->data.size() << "\n";
						//cout << "DATA" << current_node->data.size() << "\n";


						//result_number = result_tree(&new_node, 0, result_number);
					}

					if (current_node.data[i].matrix_graph[find_i][find_k] == 0)
					{
						tmp_data = decompozition(current_node.data[i], find_i, find_k);
						new_data[i] = tmp_data[0];
						new_data[i + 1] = tmp_data[1];

						//////////////////////////////////////////////////
						Node new_node = create_node(new_data);
						//new_node.count_of_trashold= current_node->count_of_trashold;

						add_Node(&new_node, &current_node);
						qu.push(new_node);
						//cout << current_number << qu.size() << "\n";
						//cout << current_node->data.size() << "\n";
						//cout << "DATA" << current_node->data.size() << "\n";


					//	result_number = result_tree(&new_node, 0, result_number);
					}

					if (current_node.data[i].matrix_graph[find_i][find_l] == 0)
					{
						tmp_data = decompozition(current_node.data[i], find_i, find_l);
						new_data[i] = tmp_data[0];
						new_data[i + 1] = tmp_data[1];

						//////////////////////////////////////////////////
						Node new_node = create_node(new_data);
						//new_node.count_of_trashold= current_node->count_of_trashold;

						add_Node(&new_node, &current_node);
						qu.push(new_node);
						//cout << current_number << qu.size() << "\n";
						//cout << current_node->data.size() << "\n";
						//result_number = result_tree(&new_node, 0, result_number);
						//cout << "DATA" << current_node->data.size() << "\n";

					}

					if (current_node.data[i].matrix_graph[find_j][find_k] == 0)
					{
						tmp_data = decompozition(current_node.data[i], find_j, find_k);
						new_data[i] = tmp_data[0];
						new_data[i + 1] = tmp_data[1];

						//////////////////////////////////////////////////
						Node new_node = create_node(new_data);
						//new_node.count_of_trashold= current_node->count_of_trashold;

						add_Node(&new_node, &current_node);
						qu.push(new_node);
						//cout << current_number << qu.size() << "\n";
						//cout << current_node->data.size() << "\n";
						//cout << "DATA" << current_node->data.size() << "\n";


						//result_number = result_tree(&new_node, 0, result_number);
					}

					if (current_node.data[i].matrix_graph[find_j][find_l] == 0)
					{
						tmp_data = decompozition(current_node.data[i], find_j, find_l);
						new_data[i] = tmp_data[0];
						new_data[i + 1] = tmp_data[1];

						//////////////////////////////////////////////////
						Node new_node = create_node(new_data);
						//new_node.count_of_trashold= current_node->count_of_trashold;

						add_Node(&new_node, &current_node);
						qu.push(new_node);
						//cout << current_number << qu.size() << "\n";
						//cout << current_node->data.size() << "\n";
						//cout << "DATA" << current_node->data.size() << "\n";


						//result_number = result_tree(&new_node, 0, result_number);
					}

					if (current_node.data[i].matrix_graph[find_k][find_l] == 0)
					{
						tmp_data = decompozition(current_node.data[i], find_k, find_l);
						new_data[i] = tmp_data[0];
						new_data[i + 1] = tmp_data[1];

						//////////////////////////////////////////////////
						Node new_node = create_node(new_data);
						//new_node.count_of_trashold= current_node->count_of_trashold;

						add_Node(&new_node, &current_node);
						qu.push(new_node);
						//cout << current_number << qu.size() << "\n";
						//cout << current_node->data.size() << "\n";
						//cout << "DATA" << current_node->data.size() << "\n";

						//result_number = result_tree(&new_node, 0, result_number);
					}
					break;

					///��������� ����� ������� � ������ � �������, count_of_trashold �������
					/// ��������� ��� �� ������� ��� ����� ������� � start_index = i
				//	cout << "DATA" << current_node->data.size() << "\n";
				}
				//cout << "DATA" << current_node->data.size() << "\n";
			}
		//	current_node.data = delete_subgraphs(&current_node);

			//if (current_node.data.size() == 4)
			//	{
			//			cout<<current_node.count_of_trashold;
				//}
			if (current_node.count_of_trashold == current_node.data.size())
				{
					
					int layer = current_node.layer;



					//->data = delete_subgraphs(current_node);
				/* if(current_node->count_of_trashold == 2)*/
					qu.pop();

					//current_node->data = delete_subgraphs(current_node);
					// cout<<"THIS";
					//current_node->data = delete_dublicate(current_node->data);
					//
					if (result > current_node.data.size())
					{
						find = 1;
						result = current_node.data.size();
						//cout << result;
						if(!(qu.empty()))
							if(qu.front().layer>layer)
								while (!qu.empty()) 
								{

									qu.pop();
								}
					

						//cout << "\nNUMBER IS CHANGE\n"<< result_number<<"\n";

						/*for (int ind_in_data = 0; ind_in_data < current_node->data.size(); ind_in_data++)
						{
							cout << "index " << ind_in_data << "\n";
							print_graph(current_node->data[ind_in_data]);
						}*/
						/*for (int ind_in_data = 0; ind_in_data < current_node->data.size(); ind_in_data++)
						{
							out_file << "index " << ind_in_data << "\n";
							current_node->data[ind_in_data].print_graph();
						}*/

						/*while (!qu.empty()) {
							qu.pop();
						}*/
					}
					else
					{
						if (!(qu.empty()))
						{	
							if((qu.front().layer>layer)&&(find == 1))
							{
								while (!qu.empty()) 
										{

											qu.pop();
										}
							}

						}

					}

			}
			else
			{
				int layer = current_node.layer;
				qu.pop();

				if (!(qu.empty()))
				{	
					if((qu.front().layer>layer)&&(find == 1))
					{
						while (!qu.empty()) 
								{

									qu.pop();
								}
					}

				}
			}
	}
	return result;
}



/*int get_threshold_number_with_decition_(graph& checked_graph)
{
	int result = 0;
	bool are_trashhold = false;
	Node root = get_root(checked_graph);
	queue<Node> qu;
	if (check_for_threshold_graph_with_subgraphs(checked_graph))
	{
		result = 1;
	}
	else
	{
		//		vector<graph> result_decompozition;
		qu.push(root);
		result = result_tree_long(&root, 0, checked_graph.size, qu);
		//		result = size(result_decompozition);
	}
	
	return result;
}*/



int result_tree(Node* current_node, int start_index_in_data, int current_number)
{
	cout << current_number << " "  << current_node->data.size() << " " << current_node->layer << "\n";

	vector<graph> result_decompozition;
	int find_i = -1, find_j=-1, find_k=-1, find_l=-1;
	/*if (current_node->root != NULL) {
		cout << "ROOT IS\n";
		for (int k = 0; k < current_node->root->data.size(); ++k)
		{
			cout << "           " << k << "\n";
			print_graph(current_node->root->data[k]);
			cout << "         \n";
		}
		cout << "\n";
	}

	cout << "DATA IS\n";


	for (int k = 0; k < current_node->data.size(); ++k)
	{
		cout << "         "<<k<<"\n";
		print_graph(current_node->data[k]);
		cout << "           \n";
	}

	cout << "NUMBER IS " << current_number<<"\n";
	cout << "CURRENT TRASHOLD NUMBER IS " << current_node->count_of_trashold << "\n";
	*/
	vector<int> point(4, -1);
	int result_number = current_number;
	for (int i = start_index_in_data; i < current_node->data.size(); ++i)
	{


		point = get_points(current_node,i);
		find_i = point[0];
		find_j = point[1];
		find_k = point[2];
		find_l = point[3];




		if (find_i==-1)
		{
			current_node->count_of_trashold++;
			//cout << "NEW COUNT IS " << current_node->count_of_trashold << "\n";
		}
		else
		{
			///////////////////////////////////////////////////////////////
			
			

			///������������, �������� ����� ������ data
			vector<graph>new_data(current_node->data.size() + 1);
			for (int j = 0; j < i; ++j)
				new_data[j] = current_node->data[j];
			for (int j = i + 1; j < current_node->data.size(); ++j)
				new_data[j + 1] = current_node->data[j];
			if (current_node->data[i].matrix_graph[find_i][find_j] == 0)
			{
				vector<graph>tmp_data = decompozition(current_node->data[i], find_i, find_j);
				new_data[i] = tmp_data[0];
				new_data[i + 1] = tmp_data[1];

				//////////////////////////////////////////////////
				Node new_node = create_node(new_data);
				//new_node.count_of_trashold= current_node->count_of_trashold;

				add_Node(&new_node, current_node);

				result_number = result_tree(&new_node, 0, result_number);
			}

			if (current_node->data[i].matrix_graph[find_i][find_k] == 0)
			{
				vector<graph>tmp_data = decompozition(current_node->data[i], find_i, find_k);
				new_data[i] = tmp_data[0];
				new_data[i + 1] = tmp_data[1];

				//////////////////////////////////////////////////
				Node new_node = create_node(new_data);
				//new_node.count_of_trashold= current_node->count_of_trashold;

				add_Node(&new_node, current_node);

				result_number = result_tree(&new_node, 0, result_number);
			}
			if (current_node->data[i].matrix_graph[find_i][find_l] == 0)
			{
				vector<graph>tmp_data = decompozition(current_node->data[i], find_i, find_l);
				new_data[i] = tmp_data[0];
				new_data[i + 1] = tmp_data[1];

				//////////////////////////////////////////////////
				Node new_node = create_node(new_data);
				//new_node.count_of_trashold= current_node->count_of_trashold;

				add_Node(&new_node, current_node);

				result_number = result_tree(&new_node, 0, result_number);
			}
			if (current_node->data[i].matrix_graph[find_j][find_k] == 0)
			{
				vector<graph>tmp_data = decompozition(current_node->data[i], find_j, find_k);
				new_data[i] = tmp_data[0];
				new_data[i + 1] = tmp_data[1];

				//////////////////////////////////////////////////
				Node new_node = create_node(new_data);
				//new_node.count_of_trashold= current_node->count_of_trashold;

				add_Node(&new_node, current_node);

				result_number = result_tree(&new_node, 0, result_number);
			}

			if (current_node->data[i].matrix_graph[find_j][find_l] == 0)
			{
				vector<graph>tmp_data = decompozition(current_node->data[i], find_j, find_l);
				new_data[i] = tmp_data[0];
				new_data[i + 1] = tmp_data[1];

				//////////////////////////////////////////////////
				Node new_node = create_node(new_data);
				//new_node.count_of_trashold= current_node->count_of_trashold;

				add_Node(&new_node, current_node);

				result_number = result_tree(&new_node, 0, result_number);
			}

			if (current_node->data[i].matrix_graph[find_k][find_l] == 0)
			{
				vector<graph>tmp_data = decompozition(current_node->data[i], find_k, find_l);
				new_data[i] = tmp_data[0];
				new_data[i + 1] = tmp_data[1];

				//////////////////////////////////////////////////
				Node new_node = create_node(new_data);
				//new_node.count_of_trashold= current_node->count_of_trashold;

				add_Node(&new_node, current_node);

				result_number = result_tree(&new_node, 0, result_number);
			}
			///��������� ����� ������� � ������ � �������, count_of_trashold �������
			/// ��������� ��� �� ������� ��� ����� ������� � start_index = i
		}
	}
	//vector<graph> new_data= delete_dublicate(current_node->data);

//cout<<"RTRW";
//vector<graph> new_data = delete_subgraphs(current_node->data);
	//current_node->data = delete_subgraphs(current_node->data);
	//current_node->data = delete_subgraphs(current_node);
	//current_node->data = delete_subgraphs(current_node);
	//current_node->data = delete_subgraphs(current_node);

	if (current_node->count_of_trashold == current_node->data.size())
	{

		
			//->data = delete_subgraphs(current_node);
		 /* if(current_node->count_of_trashold == 2)
			for (int ind_in_data = 0; ind_in_data < current_node->data.size(); ind_in_data++)
			{
				out_file << "index " << ind_in_data << "\n";
				current_node->data[ind_in_data].print_graph();
			}*/
		

		//current_node->data = delete_subgraphs(current_node);
		// cout<<"THIS";
		 //current_node->data = delete_dublicate(current_node->data);
		 //
		if (result_number > current_node->data.size())
		{
			result_number = current_node->data.size();

			//cout << "\nNUMBER IS CHANGE\n"<< result_number<<"\n";

			/*for (int ind_in_data = 0; ind_in_data < current_node->data.size(); ind_in_data++)
			{
				cout << "index " << ind_in_data << "\n";
				print_graph(current_node->data[ind_in_data]);
			}*/
		}
	}

	return result_number;
}

int get_threshold_number_with_decition_tree(graph& checked_graph)
{
	int result = 0;
	bool are_trashhold = false;
	Node root = get_root(checked_graph);
	if (check_for_threshold_graph_with_subgraphs(checked_graph))
	{
		result = 1;
	}
	else
	{
		//		vector<graph> result_decompozition;

		result = result_tree(&root, 0, checked_graph.size);
		//		result = size(result_decompozition);
	}
	return result;
}







void get_result(int size)
{
	graph gr(size);
	//gr.print_graph();
	out_file << check_for_threshold_graph_with_subgraphs(gr);

}

void get_result_for_trashold_number(int size)
{
	graph gr(size);

	//gr.print_graph();
	//out_file << get_threshold_number_with_decition_tree(gr) << " ";
	//out_file << get_threshold_number_with_decition_tree_long(gr) << " ";
	out_file << get_threshold_number_with_decition_tree_long_cycle(gr) << " \n";

}

void get_result_for_estimate_trashold_number(int size)
{
	graph gr(size);
	//gr.print_graph();
	out_file << get_estimate_threshold_number(gr)<<"\n";
}


int main()
{
	int size_of_data;
	file >> size_of_data;
	int size;
	file >> size;
	//double degree;
	//file>>degree;
	//out_file<<degree<<"\n";
	//cout << size_of_data << " " << size<<" ";
	
	for (int i = 0; i < size_of_data; ++i)
	{
		cout<<i;
		//graph gr(size);
	get_result_for_trashold_number(size);
	//get_result_for_estimate_trashold_number(size);
	}
	//get_result_for_trashold_number(size);


//system("pause");
}