#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <fstream>
#include <sstream>
#include <Eigen/Cholesky>
#include <Eigen/Sparse>
#include <string>
#include <vector>
#define REP(i,j,k) for(int i = j ; i < k ; ++i)

using namespace std;

struct Problem
{
    int mdim = 0;
    int ndim = 0;
    int nblock = 0;
    float theta = 1.0;
    vector<int> block_struct;
    vector<int> block_base;
    Eigen::VectorXf c;
    vector<Eigen::MatrixXf,Eigen::aligned_allocator<Eigen::MatrixXf> > F;

    Eigen::VectorXf xVec;
    Eigen::MatrixXf xMat;
    Eigen::MatrixXf yMat;
};

int parse_input(char* prob_path, Problem & prob) //, int & mdim, int & nblock, int &vector<int> block_struct,)
{

    ifstream prob_file(prob_path);

    if (prob_file.is_open())
    {

        string line;
        int idx = 0;
        int mdim = 0;
        int nblock = 0;


        vector<int> block_struct;
        int prev_base = 0;
        vector<int> block_base;
        while ( idx <= 2 ) 
        {
            getline (prob_file,line);
            //cout << line << '\n'; 
            if(line[0] != '"')
            {
                replace( line.begin(), line.end(), ',', ' ');
                replace( line.begin(), line.end(), '{', ' ');
                replace( line.begin(), line.end(), '}', ' ');

                istringstream iss(line);
                if (idx == 0)
                {
                    iss >> mdim;
                }
                else if (idx == 1)
                {
                    iss >> nblock;
                }
                else if (idx == 2)
                {
                    REP(i,0,nblock)
                    {
                        int block_size;
                        iss >> block_size;

                        if (block_size < 0)
                        {
                            block_size *= -1;
                        }

                        block_struct.push_back(block_size);
                        block_base.push_back(prev_base);
                        prev_base += block_size;
                    }
                }
                else 
                {
                    break;
                }
                idx += 1;
            }
        }

        //vector<Eigen::MatrixXf> F;
        Eigen::VectorXf c(mdim);
        //Eigen::VectorXf c;
        vector<Eigen::MatrixXf,Eigen::aligned_allocator<Eigen::MatrixXf> > F;


        int ndim = 0;

        REP(i,0,nblock)
        {
            ndim += block_struct[i];
        }

        //cout << ndim << endl;

        REP(i,0,mdim+1)
        {
            Eigen::MatrixXf mat(ndim,ndim);
            mat.setZero();
            F.push_back(mat);
        }


        while ( getline (prob_file,line) ) 
        { 
            //cout << line << '\n'; 
            if(line[0] != '"')
            {
                replace( line.begin(), line.end(), ',', ' ');
                replace( line.begin(), line.end(), '{', ' ');
                replace( line.begin(), line.end(), '}', ' ');

                istringstream iss(line);

                if( idx == 3)
                {
                    REP(i,0,mdim)
                    {
                        float val;
                        iss >> val;
                        c[i] = val;
                    }

                }
                else
                {
                    int mat_idx = 0;
                    int blk_idx = 0;
                    int i_idx = 0;
                    int j_idx = 0;
                    float val = 0;

                    iss >> mat_idx >> blk_idx >> i_idx >> j_idx >> val;

                    int base = block_base[blk_idx-1];

                    int i_abs = i_idx - 1 + base;
                    int j_abs = j_idx - 1 + base;

                    F[mat_idx](i_abs,j_abs) = val;
                    F[mat_idx](j_abs,i_abs) = val;

                }
                idx += 1;
            }
        }
        
        //cout << c << endl;
        //REP(i,0, mdim+1)
        //{
        //    cout << F[i] << endl;
        //}
        prob.F = F;
        prob.c = c;
        prob.mdim = mdim;
        prob.ndim = ndim;
        prob.nblock = nblock;
        prob.block_struct = block_struct;
        prob.block_base = block_base;
        prob_file.close();
    }
    {
        return -1;
    }

    return 0;
}


int parse_init(char* init_path, Problem & prob)
{
    //cout << init_path << endl;
    ifstream prob_file(init_path);

    if (prob_file.is_open())
    {
        string line;
        int nblock = prob.nblock;
        int mdim = prob.mdim;
        int ndim = prob.ndim;

        vector<int> block_struct = prob.block_struct;
        vector<int> block_base = prob.block_base;


        Eigen::VectorXf xVec(mdim);
        Eigen::MatrixXf xMat(ndim,ndim);
        Eigen::MatrixXf yMat(ndim,ndim);
        xVec.setZero();
        xMat.setZero();
        yMat.setZero();

        int idx = 0;
        
        while ( getline (prob_file,line) ) 
        {

            replace( line.begin(), line.end(), ',', ' ');
            replace( line.begin(), line.end(), '{', ' ');
            replace( line.begin(), line.end(), '}', ' ');
            istringstream iss(line);

            if (idx == 1)
            {
                //parse xVec
                float val;
                REP(i, 0, mdim)
                {
                    iss >> val;
                    xVec[i] = val;
                }

                //cout << xVec << endl;
            }
            else if ( idx >=3 && idx <= 2 + nblock )
            {
                //parse xMat
                int blk_nth = idx - 3;
                float val;
                int size = block_struct[blk_nth];
                int base = block_base[blk_nth];

                REP(i, 0, size)
                {
                    REP(j, 0, size) 
                    { 
                        iss >> val; 
                        int abs_i = i + base;
                        int abs_j = j + base;
                        xMat(abs_i,abs_j) = val;
                    }
                }
            }
            else if (idx >= 2 + nblock + 2)
            {
                //parse yMat
                int blk_nth = idx - (2+nblock+2);
                float val;
                int size = block_struct[blk_nth];
                int base = block_base[blk_nth];

                REP(i, 0, size)
                {
                    REP(j, 0, size) 
                    { 
                        iss >> val; 
                        int abs_i = i + base;
                        int abs_j = j + base;
                        yMat(abs_i,abs_j) = val;
                    }
                }
            }

            idx += 1;
        }
        prob.xVec=xVec;
        prob.xMat=xMat;
        prob.yMat=yMat;

        //cout << xVec << endl;
        //cout << xMat << endl;
        //cout << yMat << endl;

    }
    else
    {
        return 1;
    }
    

    return 0;
}




bool is_psd(Eigen::MatrixXf X)
{
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es(X);
    Eigen::VectorXf eivals = es.eigenvalues();
    int len = eivals.size();

    REP(i, 0, len)
    {
        if (eivals(i) <= 0)
        {
            return false;
        }
    }

    return true;
}

Eigen::VectorXf flatten(Eigen::MatrixXf mat)
{
    Eigen::MatrixXf temp = mat;
    temp.transposeInPlace();

    Eigen::VectorXf ret(Eigen::Map<Eigen::VectorXf>(temp.data(), temp.cols()*temp.rows()));
    return ret;
}



Eigen::MatrixXf concatenate( vector<Eigen::MatrixXf,Eigen::aligned_allocator<Eigen::MatrixXf> > & list, int axis = 0 )
{
    int n = list.size();
    int rows = 0;
    int cols = 0;

    int each_rows = 0;
    int each_cols = 0;
    
    vector<int> dims;
    if( axis == 0 )
    {
        REP(i,0,n)
        {
            each_rows = list[i].rows();
            each_cols = list[i].cols();

            rows += each_rows;
            dims.push_back(each_rows);
            cols = each_cols;
        }

        Eigen::MatrixXf ret(rows, cols);
        ret.setZero();

        int start_dim = 0;
        REP(i,0,n) 
        { 
            int i0 = start_dim;
            int j0 = 0;
            ret.block(i0,j0,dims[i],cols) = list[i];
            start_dim += dims[i];
        }

        return ret;
    }
    else
    {
        REP(i,0,n)
        {
            each_rows = list[i].rows();
            each_cols = list[i].cols();

            rows = each_rows;
            cols += each_cols;
            dims.push_back(each_cols);
        }

        Eigen::MatrixXf ret(rows, cols);
        ret.setZero();

        int start_dim = 0;
        REP(i,0,n) 
        { 
            int i0 = 0;
            int j0 = start_dim;
            ret.block(i0,j0,rows,dims[i]) = list[i];
            start_dim += dims[i];
        }

        return ret;
    }
}

Eigen::MatrixXf linear_prod( Eigen::VectorXf & y, vector<Eigen::MatrixXf,Eigen::aligned_allocator<Eigen::MatrixXf> > & A)
{
    Eigen::MatrixXf ret = A[0];
    int n = A.size();

    REP(i,1,n)
    {
        ret += y(i) * A[i];
    }

    return ret;
}


// Reference paper: http://ocw.mit.edu/courses/electrical-engineering-and-computer-science/6-251j-introduction-to-mathematical-programming-fall-2009/readings/MIT6_251JF09_SDP.pdf
//
// Solve the normal equations for D and y:
// C - theta*X^-1 + theta * X^-1 * D * X^-1 = sum_{i=1}^m y_i A_i
// A_i*D = 0, i = 1,...,m
float solve_normal_equation(
        int mdim,
        int ndim,
        Eigen::MatrixXf &X, 
        vector<Eigen::MatrixXf,Eigen::aligned_allocator<Eigen::MatrixXf> > & A,
        float theta,
        Eigen::MatrixXf &C,
        Eigen::MatrixXf &D,
        Eigen::VectorXf &y)
{



    vector<Eigen::MatrixXf,Eigen::aligned_allocator<Eigen::MatrixXf> > temp;
    temp.push_back(flatten(X-1/theta*X*C*X));

    Eigen::MatrixXf temp_zeros(mdim,1);
    temp_zeros.setZero();
    temp.push_back(temp_zeros);
    Eigen::MatrixXf N = concatenate(temp);
    //cout << N << endl;
    temp.clear();


    REP(i,0,mdim)
    {
       //cout << "before" << endl;
       //cout <<  -1/theta *X*A[i]*X << endl;
       //cout << "after" << endl;
       //cout <<  flatten(-1/theta *X*A[i]*X) << endl;
       temp.push_back(flatten(-1/theta *X*A[i]*X));
    }


    Eigen::MatrixXf Q = concatenate(temp,1);

    
    //cout << Q << endl;
    temp.clear();

    REP(i,0,mdim)
    {
       //cout << "before" << endl;
       //cout << A[i] << endl;
       //cout << "after" << endl;
       //cout << flatten(A[i]).tranpse() << endl;
       temp.push_back(flatten(A[i]).transpose());
    }

    //cout << A[mdim-2] << endl;
    Eigen::MatrixXf P = concatenate(temp);
    temp.clear();


    Eigen::MatrixXf M(ndim*ndim + mdim, ndim*ndim + mdim);
    M.setZero();

    Eigen::MatrixXf ind(ndim*ndim,ndim*ndim);
    Eigen::MatrixXf zeros(mdim,mdim);
    ind.setIdentity();
    zeros.setZero();
    M.block(0,0,ndim*ndim,ndim*ndim) = ind;
    M.block(0,ndim*ndim, ndim*ndim, mdim) = Q;
    M.block(ndim*ndim,0, mdim, ndim*ndim) = P;
    M.block(ndim*ndim,ndim*ndim, mdim,mdim) = zeros;


    
    //cout << M << endl;
    Eigen::VectorXf Z = M.colPivHouseholderQr().solve(N);

    float err = (M*Z-N).norm();



    D = Z.block(0,0,ndim*ndim,1);
    y = Z.block(ndim*ndim,0, mdim,1);


    D.resize(ndim,ndim);
    D = 0.5 * (D + D.transpose());
    
    //cout << M << endl;
    //cout << N << endl;
    //cout << D << endl;
    //cout << y << endl;

    return err;

}


float prod_sum(Eigen::MatrixXf & x, Eigen::MatrixXf &y)
{

    float ret = 0.0;

    int rows = x.rows();
    int cols = x.cols();


    REP(i,0,rows)
    {
        REP(j,0,cols)
        {
            ret += x(i,j) * y(i,j);
        }
    }

    return ret;
}


void sdp(Problem & prob)
{
    int mdim = prob.mdim;
    int ndim = prob.ndim;
    int nblock = prob.nblock;
    float theta = prob.theta;
    vector<int> & block_struct = prob.block_struct;
    vector<int> & block_base = prob.block_base;

    Eigen::MatrixXf X = prob.yMat;
    Eigen::MatrixXf C = -prob.F[0];
    Eigen::VectorXf y;
    Eigen::MatrixXf S;
    Eigen::VectorXf b = -prob.c;

    vector<Eigen::MatrixXf,Eigen::aligned_allocator<Eigen::MatrixXf> > A;

    REP(i,0,mdim)
    {
        Eigen::MatrixXf mat = -prob.F[i+1];
        A.push_back(mat);
    }


    int m = mdim;
    int n = ndim;

    while (true)
    {

        float alpha = 0.8;
        theta = alpha * theta;

        Eigen::MatrixXf D;
        float err= solve_normal_equation( mdim, ndim, X, A, theta, C, D, y);

        if (err > 1e-2)
        {
            break;
        }

        S = C - linear_prod(y,A);


        float t = 1.0;

        while( !is_psd(X + t*D))
        {
            t = t * alpha;
        }

        float primal_criterion = -b.transpose()*y;
        float dual_criterion = -prod_sum(C,X);

        cout << "primal_criterion: " << primal_criterion << endl;
        cout << "dual_criterion: " << dual_criterion << endl;

        if (theta < 1e-4)
        { 
            break;
        }
        else
        { 
            X = X + t*D;
        }

    }

    prob.yMat = X;
    prob.xMat = S;
    prob.xVec = y;
}


void output_sol(Problem & prob, char * path)
{
    ofstream output(path);

    output << "xVec" << endl;
    output << prob.xVec << endl;
    output << "xMat" << endl;
    output << prob.xMat << endl;
    output << "yMat" << endl;
    output << prob.yMat << endl;
    output.close();
}

int main(int argc, char *argv[])
{
    if (argc != 4)
    {
        cout << "Usage: sdp problem_file init_file output_file" << endl;
        exit(1);
    }

    char * prob_path = argv[1];
    char * init_path = argv[2];

    Problem prob;
    cout << "Parsing problem..." << endl;
    parse_input(prob_path, prob);
    cout << "Parsing initialization..." << endl;
    parse_init(init_path, prob);

    cout << "Solving..." << endl;
    sdp(prob);
    cout << "Writing down solution..." << endl;
    output_sol(prob, argv[3]);
}
