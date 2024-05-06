#include "PSToolboxRunner.h"
#include "PSToolboxBaseEdge.h"
#include "Connector.h"
#include "my_tools.h"
#include <Eigen/Dense>
#include <fstream>
#include <iostream>

static vector<vector <int>> DEFAULT_VECTOR;

PSToolboxRunner::PSToolboxRunner(
    vector<PSToolboxBaseEdge *> &_e, 
    vector<Connector *> &_c,
    vector<int> &_con_at_edge_start,
    vector<int> &_con_at_edge_end,
	vector< vector <int> > & _rigid_subsystems = DEFAULT_VECTOR){
  edges = _e;
  cons = _c;
  con_at_edge_start = _con_at_edge_start;
  con_at_edge_end = _con_at_edge_end;
  save_data=false;
  folder = "data";
  DEBUG=false;
  Node_mul=1;
  set_limit = false;
  p_limMin = 0.0;
  p_limMax = 0.0;
  write_interval = 1000.0;
  save_interval = 1.0;
  p_start_time = 0.0;
  p_end_time = 1.0e10;
  for (int i=0; i<_rigid_subsystems.size(); i++)
      rigid_subsystems.push_back(_rigid_subsystems.at(i));

	// Build info of rigid subsystems

	int idx;
	vector<int> tmp;
	vector<bool> tmp_bool;
	vector< vector<int> > tmp2;
	vector< vector<bool> > tmp_bool2;
	vector <int> edge_list, edge_tmp;

	DEBUG=false;

	for (int i=0; i<rigid_subsystems.size(); i++){
		edge_list.clear();
		tmp2.clear();
		tmp_bool2.clear();
		edge_tmp.clear();

		// Collect node indices
		tmp.clear();
		for (int j=0; j<rigid_subsystems.at(i).size(); j++){
			idx = rigid_subsystems.at(i).at(j);
			if (con_at_edge_start.at(idx)>-1)
				tmp.push_back(con_at_edge_start.at(idx));
			if (con_at_edge_end.at(idx)>-1)
				tmp.push_back(con_at_edge_end.at(idx));
		}
		sort( tmp.begin(), tmp.end() );
		tmp.erase( unique( tmp.begin(), tmp.end() ), tmp.end() );
		rs_idx_nodes.push_back(tmp);


		// Now build connectivity info
		for (int j=0; j<rs_idx_nodes.at(i).size(); j++){
			tmp.clear();
			tmp_bool.clear();
			idx = rs_idx_nodes.at(i).at(j);

			for (int k=0; k<edges.size(); k++){
				if (con_at_edge_end.at(k)==idx){
					tmp.push_back(k);
					tmp_bool.push_back(false);
					edge_list.push_back(k);
				}
			}

			for (int k=0; k<edges.size(); k++){
				if (con_at_edge_start.at(k)==idx){
					//cout<<endl<<"edges.at(i)->name: "<<edges.at(k)->name;
					if (!(edges.at(k)->edge_type=="Tank")){
						tmp.push_back(k);
						tmp_bool.push_back(true);
						edge_list.push_back(k);
					}
				}
			}

			tmp2.push_back(tmp);
			tmp_bool2.push_back(tmp_bool);
			edge_tmp.insert(end(edge_tmp), begin(tmp), end(tmp));
		}
		rs_idx_nodes_edges.push_back(tmp2);
		rs_idx_nodes_is_front.push_back(tmp_bool2);

		sort( edge_list.begin(), edge_list.end() );
		edge_list.erase( unique( edge_list.begin(), edge_list.end() ), edge_list.end() );
		num_of_edges.push_back(edge_list.size());

		sort( edge_tmp.begin(), edge_tmp.end() );
		edge_tmp.erase( unique( edge_tmp.begin(), edge_tmp.end() ), edge_tmp.end() );
		rs_idx_edges.push_back(edge_tmp);

	}

	for (int i=0; i<rigid_subsystems.size(); i++){

		if (DEBUG){
			cout<<endl<<"Topology of of subsystem "<<i;
			cout<<endl<<"\t # of nodes: "<<rs_idx_nodes.at(i).size();
			cout<<", # of edges:"<<num_of_edges.at(i);
			cout<<endl<<"\t list of nodes: ";
			for (auto ii: rs_idx_nodes.at(i))
				std::cout << ii << " ";

			cout<<endl<<"\t list of edges: ";
			for (auto ii: rs_idx_edges.at(i))
				std::cout << ii << " ";

			int edge_idx;
			for (int jj=0; jj<rs_idx_nodes.at(i).size(); jj++){
				edge_idx=rs_idx_nodes.at(i).at(jj);
				cout << endl<<"\t\t node: "<<edge_idx<<", edges: ";
				for (int kk=0; kk<rs_idx_nodes_edges.at(i).at(jj).size(); kk++){
					cout << " "<<rs_idx_nodes_edges.at(i).at(jj).at(kk);
					if (rs_idx_nodes_is_front.at(i).at(jj).at(kk))
						cout<<"+";
					else
						cout<<"-";
				}
			}
			cout<<endl;
			cin.get();
		}
	}

}

void PSToolboxRunner::RS_Solve(int rs_idx, double t_target, bool rs_save)
{
	int n_n=rs_idx_nodes.at(rs_idx).size();
	int n_e=rs_idx_edges.at(rs_idx).size();
	MatrixXd A = MatrixXd::Zero(n_n+n_e, n_n+n_e);
	VectorXd b = VectorXd::Zero(n_n+n_e);
	VectorXd x = VectorXd::Zero(n_n+n_e);
	VectorXd xnew = VectorXd::Zero(n_n+n_e);

	// Set initial condition for faster convergence
	for (int i_edge=0; i_edge<n_e; i_edge++)
	{
		int idx_edge = rs_idx_edges.at(rs_idx).at(i_edge);
		if (edges.at(i_edge)->is_rigid_element){
			x(i_edge)=3600*edges.at(idx_edge)->Get_Q();
		}
	}
	for (int i_edge=0; i_edge<n_e; i_edge++)
	{
		int idx_edge = rs_idx_edges.at(rs_idx).at(i_edge);
		int idx_p1_rs = getIndexOf(rs_idx_nodes.at(rs_idx), con_at_edge_start.at(idx_edge));
		int idx_p2_rs = getIndexOf(rs_idx_nodes.at(rs_idx), con_at_edge_end.at(idx_edge));

		if (edges.at(idx_edge)->is_rigid_element){
			if (!(edges.at(idx_edge)->edge_type=="Tank"))
				x(n_e+idx_p1_rs)=(edges.at(idx_edge)->Get_p1())/1000./9.81;

			x(n_e+idx_p2_rs)=(edges.at(idx_edge)->Get_p2())/1000./9.81;
		}
	}
	xnew=x;
	// initialization done.

	double x_norm, xnew_norm, err=1.e5, err_max=0.001;
	int iter=0, iter_max=2000;

	if (DEBUG)
		cout<<endl<<endl<<"Updating rigid subsystem #"<<rs_idx;

	while (err>err_max)
	{
		x_norm=x.squaredNorm();
		if (iter==iter_max){
			cout<<endl<<endl<<"ERROR!!!!!!!!!!!!";
			cout<<endl<<"PSToolboxRunner::Solve() iter=iter_max!!!!";
			cout<<endl<<"Last step";

			printf("\n\t iter #%2d: x_norm=%5.2e, xnew_norm=%5.2e, err=%5.3e, err_max=%5.2e\n",iter,x_norm,xnew_norm,err,err_max);
			for (int i_edge=0; i_edge<n_e; i_edge++){
				int idx_edge = rs_idx_edges.at(rs_idx).at(i_edge);
				int idx_p1 = getIndexOf(rs_idx_nodes.at(rs_idx), con_at_edge_start.at(idx_edge));
				int idx_p2 = getIndexOf(rs_idx_nodes.at(rs_idx), con_at_edge_end.at(idx_edge));
				cout<<endl<<"\tEdge "<<rs_idx_edges.at(rs_idx).at(i_edge);
				if (idx_p1<0)
					printf(": p1=------ mwc, p2=%5.2f mwc, Q=%5.2f m3/h",x(n_e+idx_p2),x(i_edge));
				else{
					if (idx_p2<0)
						printf(": p1=%5.2f mwc, p2=------ mwc, Q=%5.2f m3/h",x(n_e+idx_p1),x(i_edge));
					else
						printf(": p1=%5.2f mwc, p2=%5.3f mwc, Q=%5.2f m3/h",x(n_e+idx_p1),x(n_e+idx_p2),x(i_edge));
				}
			}

			exit(-1);
		}
		for (int i_node=0; i_node<n_n; i_node++){
			for (int i_edge=0; i_edge<rs_idx_nodes_edges.at(rs_idx).at(i_node).size(); i_edge++){
				int idx_of_edge_at_node=rs_idx_nodes_edges.at(rs_idx).at(i_node).at(i_edge);
				int idx_of_Q_tmp=getIndexOf(rs_idx_edges.at(rs_idx), idx_of_edge_at_node);
				//idx_of_Q.push_back(idx_of_Q_tmp);

				if (idx_of_Q_tmp<0){
					cout<<endl<<"?????????????? PSToolboxRunner::Solve ??? idx_of_Q="<<idx_of_Q_tmp<<endl;
					exit(-1);
				}

				// Continuity
				if (rs_idx_nodes_is_front.at(rs_idx).at(i_node).at(i_edge))
					A(i_node,idx_of_Q_tmp)=1./3600;
				else
					A(i_node,idx_of_Q_tmp)=-1./3600;
			}
			int con_idx=rs_idx_nodes.at(rs_idx).at(i_node);
			double tmp = cons.at(con_idx)->demand;
			b(i_node)=tmp/3600;
			//cout<<endl<<"  connector #"<<con_idx<<", demand:"<<tmp;
			//cin.get();
		}

		int idx_p1, idx_p2, idx_edge;
		double LHS, coeff_Q, coeff_p1, coeff_p2;
		bool is_front_l, is_back_l;
		for (int i_edge=0; i_edge<n_e; i_edge++){
			idx_edge = rs_idx_edges.at(rs_idx).at(i_edge);
			idx_p1 = getIndexOf(rs_idx_nodes.at(rs_idx), con_at_edge_start.at(idx_edge));
			idx_p2 = getIndexOf(rs_idx_nodes.at(rs_idx), con_at_edge_end.at(idx_edge));

			is_front_l=true; is_back_l=true;
			if (idx_p1<0)
				is_front_l=false;
			if (idx_p2<0)
				is_back_l=false;

			if (!(edges.at(i_edge)->is_rigid_element))
				if (!(is_front_l || is_back_l)){
					cout<<endl<<"PSToolboxRunner::Solve() -> WTF????"<<endl;
					exit(-1);
				}

			edges.at(idx_edge)->GetEdgeEquationCoeffs(t_target, is_front_l, LHS, coeff_Q, coeff_p1, coeff_p2);

			A(n_n+i_edge,i_edge)       = coeff_Q/1000/9.81/3600;
			if (idx_p1>-1)
				A(n_n+i_edge,n_e+idx_p1) = coeff_p1;
			if (idx_p2>-1)
				A(n_n+i_edge,n_e+idx_p2) = coeff_p2;
			b(n_n+i_edge)              = LHS/1000/9.81;

			//			if (DEBUG){
			//				cout<<endl<<" edge:"<<idx_edge<<": "<<idx_p1<<" -> "<<idx_p2;
			//				if (idx_p1>-1)
			//					cout<<", p"<<rs_idx_nodes.at(rs_idx).at(idx_p1);
			//				else
			//					cout<<",  x ";
			//				if (idx_p2>-1)
			//					cout<<" -> p"<<rs_idx_nodes.at(rs_idx).at(idx_p2);
			//				else
			//					cout<<" -> x";
			//			}
		}

		//std::cout << "Here is the matrix A:\n" << A << std::endl;
		//std::cout << "Here is the vector b:\n" << b << std::endl;
		xnew = A.colPivHouseholderQr().solve(b);
		//xnew = A.householderQr().solve(b);
		double RELAX=0.85;
		xnew=RELAX*x+(1-RELAX)*xnew;
		//std::cout << "The solution is:\n" << x << std::endl;
		xnew_norm=xnew.squaredNorm();
		err=fabs(xnew_norm-x_norm)/xnew_norm;


		if (DEBUG){
			printf("\n\t iter #%2d: x_norm=%5.2e, xnew_norm=%5.2e, err=%5.3e, err_max=%5.2e\n",iter,x_norm,xnew_norm,err,err_max);
			for (int i_edge=0; i_edge<n_e; i_edge++){
				int idx_edge = rs_idx_edges.at(rs_idx).at(i_edge);
				int idx_p1 = getIndexOf(rs_idx_nodes.at(rs_idx), con_at_edge_start.at(idx_edge));
				int idx_p2 = getIndexOf(rs_idx_nodes.at(rs_idx), con_at_edge_end.at(idx_edge));
				cout<<endl<<"\tEdge "<<rs_idx_edges.at(rs_idx).at(i_edge);
				if (idx_p1<0)
					printf(": p1=------ mwc, p2=%5.2f mwc, Q=%5.2f m3/h",x(n_e+idx_p2),x(i_edge));
				else{
					if (idx_p2<0)
						printf(": p1=%5.2f mwc, p2=------ mwc, Q=%5.2f m3/h",x(n_e+idx_p1),x(i_edge));
					else
						printf(": p1=%5.2f mwc, p2=%5.3f mwc, Q=%5.2f m3/h",x(n_e+idx_p1),x(n_e+idx_p2),x(i_edge));
				}
			}
		}

		//		cout<<endl<<endl<<"\t list of nodes: ";
		//			for (auto ii: rs_idx_nodes.at(rs_idx))
		//				std::cout << ii << " ";
		//
		//			cout<<endl<<"\t list of edges: ";
		//			for (auto ii: rs_idx_edges.at(rs_idx))
		//				std::cout << ii << " ";
		//
		for (int i_edge=0; i_edge<n_e; i_edge++){
			idx_edge = rs_idx_edges.at(rs_idx).at(i_edge);
			if (edges.at(i_edge)->is_rigid_element){
				edges.at(idx_edge)->Set_Q(x(i_edge)/3600.);
			}
		}


		// for the rigid elements, pressure must be stored element-wise
		for (int i_edge=0; i_edge<n_e; i_edge++){
			idx_edge = rs_idx_edges.at(rs_idx).at(i_edge);
			int idx_p1_rs = getIndexOf(rs_idx_nodes.at(rs_idx), con_at_edge_start.at(idx_edge));
			int idx_p2_rs = getIndexOf(rs_idx_nodes.at(rs_idx), con_at_edge_end.at(idx_edge));
			//cout<<endl<<"i_edge="<<i_edge<<", idx_edge="<<idx_edge;
			//cout<<", is_rigid_element:"<<edges.at(idx_edge)->is_rigid_element;

			if (edges.at(idx_edge)->is_rigid_element){
				//cout<<", idx_p1="<<con_at_edge_start.at(idx_edge)<<" -> "<<idx_p1_rs<<"th";
				//cout<<", idx_p2="<<con_at_edge_end.at(idx_edge)<<" -> "<<idx_p2_rs<<"th";

				if (!(edges.at(idx_edge)->edge_type=="Tank"))
					edges.at(idx_edge)->Set_p1(x(n_e+idx_p1_rs)*1000.*9.81);

				edges.at(idx_edge)->Set_p2(x(n_e+idx_p2_rs)*1000.*9.81);
			}
		}
		//cin.get();

		x=xnew;
		iter++;

		//if (DEBUG){
		//	printf("\n\t iter #%2d: x_norm=%5.2e, xnew_norm=%5.2e, err=%5.3e, err_max=%5.2e",iter,x_norm,xnew_norm,err,err_max);
		//cin.get();
	}
	//cout << "iter : " << iter << endl;

	if (rs_save){
		for (int i_edge=0; i_edge<n_e; i_edge++){
			int idx_edge = rs_idx_edges.at(rs_idx).at(i_edge);
			edges.at(idx_edge)->Save_data();
		}
	}

//	for (int i_edge=0; i_edge<n_e; i_edge++){
//		int idx_edge = rs_idx_edges.at(rs_idx).at(i_edge);
//		edges.at(idx_edge)->Set_t(t_target);
//	}

	if (DEBUG){
		printf("\n\t iter #%2d: x_norm=%5.2e, xnew_norm=%5.2e, err=%5.3e, err_max=%5.2e\n",iter,x_norm,xnew_norm,err,err_max);
		for (int i_edge=0; i_edge<n_e; i_edge++){
			int idx_edge = rs_idx_edges.at(rs_idx).at(i_edge);
			int idx_p1 = getIndexOf(rs_idx_nodes.at(rs_idx), con_at_edge_start.at(idx_edge));
			int idx_p2 = getIndexOf(rs_idx_nodes.at(rs_idx), con_at_edge_end.at(idx_edge));

			cout<<endl<<"\tEdge "<<idx_edge<<" ("<<edges.at(idx_edge)->name<<")";

			if (idx_p1<0)
				printf(": p1=------ mwc, p2=%5.2f mwc, Q=%5.2f m3/h",x(n_e+idx_p2),x(i_edge));
			else{
				if (idx_p2<0)
					printf(": p1=%5.2f mwc, p2=------ mwc, Q=%5.2f m3/h",x(n_e+idx_p1),x(i_edge));
				else
					printf(": p1=%5.2f mwc, p2=%5.3f mwc, Q=%5.2f m3/h",x(n_e+idx_p1),x(n_e+idx_p2),x(i_edge));
			}
		}
		//cin.get();
	}

}

void PSToolboxRunner::Run(double t_max){

  //create vectors for critical values
  vector<double> p_crit;
  vector<double> x_crit;
  vector<double> t_crit;
  vector<double> last_save;
  vector<double> last_save_rigid;
  vector<string> ID_crit;

  // Initialization
  for (unsigned int i=0; i<edges.size(); i++)
  {
    edges.at(i)->Ini(Node_mul);
    last_save.push_back(0.0);
  }
  for (unsigned int i=0; i<rigid_subsystems.size(); i++)
  {
    last_save_rigid.push_back(0.0);
  }
  int rigid_elem_update = true;

  // Simulation
  double t_global=0., dt_out=t_max/100., t_out=-1.e-10;
  double t_next;
  double t_save = write_interval;
  int update_idx;
  vector<bool> update_edges(edges.size());
  int STEP=0;

  double T_TOL=1.e-10;
  /*for (unsigned int i=0; i<edges.size(); i++)
    if (edges.at(i)->Get_dt()<T_TOL)
      T_TOL=edges.at(i)->Get_dt();
  T_TOL/=100.;*/

  double last_PRV_step = 0.0;

  int pump_state = 0;
  double pump_last_t = 0;

  while (t_global<t_max)
  {
    /*
    ======================== APPLY RUNTIME MODIFIERS ===============================
    */

   for(int i = 0; i < modifiers.size(); i++)
   {
    if(modifiers[i].done == false && t_global >= modifiers[i].time)
    {
       if(!modifiers[i].has_started)
	   {
		cout << "\n--> Modifier started at t = " << t_global << "s!\n";
		//get old demands
		for (int j = 0; j < modifiers[i].connectors.size(); j++)
		{
			int idx = modifiers[i].connectors[j];
			modifiers[i].old_demands.push_back(cons[idx]->demand);
		}
		modifiers[i].has_started = true;
	   }

      

      bool autoset_final = false;
	  if(t_global >= modifiers[i].time + modifiers[i].set_time)
	  {
		cout << "\n--> Modifier finalized at t = " << t_global << "s!\n";
		modifiers[i].done = true;
		autoset_final = true;
	  }


      for (int j = 0; j < modifiers[i].connectors.size(); j++)
      {
        int idx = modifiers[i].connectors[j];
        double finaldemand = modifiers[i].demands[j];
		double olddemand = modifiers[i].old_demands[j];
		double m = (finaldemand - olddemand) / modifiers[i].set_time; 
		double newdemand = olddemand + m * (t_global - modifiers[i].time);

		if(autoset_final) newdemand = finaldemand;

        //connecting multiple pipes
        if(cons[idx]->type == 0)
        {
          cons[idx]->demand = newdemand;
          //cout << " " << cons[idx]->name << " demand " << olddemand << "-> " << newdemand  << endl;  
        }      

        //connecting velocity -- pipe
        if(cons[idx]->type == 1 && cons[idx]->BC_type == "Velocity")
        {
          cons[idx]->demand = newdemand;
          //cout << " " << cons[idx]->name << " demand " << olddemand << "-> " << newdemand  << endl;  

          double vold = cons[idx]->BC_value;
          double vnew = newdemand / olddemand * vold;
          cons[idx]->BC_value = vnew;
          //cout << " " << cons[idx]->name << " velocity " << vold << "-> " << vnew  << endl;
        }

        //connecting pressure -- pipe
        if(cons[idx]->type == 1 && cons[idx]->BC_type == "Pressure")
        {
          cons[idx]->BC_value += newdemand;
          //cout << " " << cons[idx]->name << " pressure -> " << cons[idx]->BC_value  << endl;
		  if(modifiers[i].set_time > 0)
		  {
			cout << "set_time cannot be > 0 for pressure BC-s! NOT IMPEMENTED";
			while(1) 1;
		  }
        }
      }
    }
   }

    /*
    ======================== TIMESTEP PREPARATIONS ===============================
    */
    STEP++;
    if (DEBUG)
    {
      cout<<endl<<"STEP     : "<<STEP;
      cout<<endl<<"t_global : "<<t_global;
    }
    if (t_global>t_out)
    {
      cout<<endl<<round(t_global/t_max*100)<<"%";
      t_out+=dt_out;
    }

    // Find smallest next target time
    fill(update_edges.begin(),update_edges.end(),false);

    update_idx=0;
    t_next=1.e5;

      if (DEBUG)
      cout<<endl<<"Timestep selection (flexible elements only)...";
		  for (unsigned int i=0; i<edges.size(); i++){
			if (DEBUG){
				if (!edges.at(i)->is_rigid_element)
					printf("\n\t edge #%2d %10s: t=%5.3e, tnext=%5.3e, current tnext=%5.3e",
							i,edges.at(i)->Get_name().c_str(),edges.at(i)->Get_t(),edges.at(i)->Get_tnext(),t_next);
      }
      if (!edges.at(i)->is_rigid_element)
      if (edges.at(i)->Get_tnext()<t_next)
      {
        update_idx=i;
        t_next=edges.at(i)->Get_tnext();

        if (DEBUG)
          cout<<"  <-- t_target ";
      }
    }
    update_edges.at(update_idx)=true;

    // Now find all elements within T_TOL
    for (unsigned int i=0; i<edges.size(); i++){
			if (!edges.at(i)->is_rigid_element)
      if (fabs((edges.at(i)->Get_tnext())-t_next)<T_TOL)
        update_edges.at(i)=true;
    }
    if (DEBUG)
    {
      cout<<endl<<endl<<"The following edges will be updated: ";
      for (unsigned int i=0; i<edges.size(); i++)
      {
        if (update_edges.at(i))
          cout<<" "<<i;
      }
    	cin.get();
	}

    /*
	======================== UPDATE RIGID SUBSYSTEMS ===============================
	*/

	if(rigid_elem_update)
	{
		for (int irs=0; irs<rigid_subsystems.size(); irs++)
		{
			bool rs_save_data;
			if(last_save_rigid[irs] + save_interval < t_next)
			{
				rs_save_data = true;
				last_save_rigid[irs] = int(t_next/save_interval)*save_interval;
			}
			else
			{
				rs_save_data = false;
			}
			Update_RS(irs,t_next,rs_save_data);
		}
		if (DEBUG)
		{
			cin.get();
		}
	}
	rigid_elem_update = false;

    /*
    ======================== UPDATE EDGES ===============================
    */
    
    for (unsigned int i=0; i<edges.size(); i++)
    { 
      if (update_edges.at(i))
      {
        if (DEBUG)
          cout<<endl<<"\t Updating edge "<<i<<" ... ";
        
		if (DEBUG)
			cout<<endl<<"\t\t UpdateInternal()... ";

		
		edges.at(i)->UpdateInternal(t_next);
		
		cons.at(con_at_edge_start.at(i))->Update(t_next,i);
		cons.at(con_at_edge_end.at(i))->Update(t_next,i);
	

        if (DEBUG)
        {  
          cout<<"done.";
          //cout << edges.at(i)->Info();
        }

        edges.at(i)->UpdateTime(t_next);

        if (save_data)
        {
          if(last_save[i] + save_interval < t_next)
          {
            edges.at(i)->Save_data();
            last_save[i] = int(t_next/save_interval)*save_interval;
          }
        }
          
        if(set_limit && t_next > p_start_time && t_next < p_end_time)
        {
          edges.at(i)->GetLargePressureValues(p_limMax,x_crit,p_crit,t_crit,ID_crit);
          edges.at(i)->GetSmallPressureValues(p_limMin,x_crit,p_crit,t_crit,ID_crit);
        }

		//rigid element update?
		auto con_start = cons.at(con_at_edge_start.at(i));
		auto con_end = cons.at(con_at_edge_end.at(i));
		if(con_start -> connected_to_rigid())
		{
			//cout << "Rigid connection, node " << con_start->name << endl;
			rigid_elem_update = true;
		}
		if(con_end -> connected_to_rigid())
		{
			//cout << "Rigid connection, node " << con_end->name << endl;
			rigid_elem_update = true;
		}

      }//end of if - update edges
    }//end of for - update edges

    if (DEBUG)
    	cin.get();

		/*
		======================== ADVANCE TIME ===============================
		 */

    t_global=t_next;

    /*
    ======================== Save at end of timestep ===============================
    */

    if (save_data && t_global > t_save)
    {
      t_save = t_next + write_interval;
      for (unsigned int i=0; i<edges.size(); i++)
      {
          edges.at(i)->Write_data(folder);
      }
    }


  }//end of while 


  /*
  ======================== Save at end of simulation ===============================
  */
  if (save_data)
  {
    for (unsigned int i=0; i<edges.size(); i++)
    {
        edges.at(i)->Write_data(folder);
    }
  }
      

  if(set_limit)
  {
    ofstream ofsMax = ofstream(limits_fname + "_max.log");
    ofstream ofsMin = ofstream(limits_fname + "_min.log");
    if(limits_type == "All")
    {
      for (int i = 0; i < ID_crit.size(); i++)
      {
        if (p_limMax < p_crit.at(i))
        {
          ofsMax << "Pressure limit exceeded in pipe " << ID_crit.at(i);
          ofsMax << "\t t = "<< t_crit.at(i) << " s";
          ofsMax << "\t x = "<< x_crit.at(i) << " m";
          ofsMax << "\t p = "<< 1e-5*p_crit.at(i) << " bar" << endl;
        }
        else
        {
          ofsMin << "Pressure limit exceeded in pipe " << ID_crit.at(i);
          ofsMin << "\t t = "<< t_crit.at(i) << " s";
          ofsMin << "\t x = "<< x_crit.at(i) << " m";
          ofsMin << "\t p = "<< 1e-5*p_crit.at(i) << " bar" << endl;
        }
      } 
    }
    if(limits_type == "Maximum")
    {
      for(int i = 0; i < edges.size(); i++) //go through all the pipes
      {
        int countMax = 0;
        int countMin = 0;
        double p_max = p_limMax;
        double p_min = p_limMin;
        int idxMax = -1;
        int idxMin = -1;
        for (int j = 0; j < ID_crit.size(); j++) //go through dataset
        {
          if(ID_crit.at(j) == edges.at(i)->Get_name() && p_crit.at(j) >= p_max)
          {
            p_max = p_crit.at(j);
            idxMax = j;
            countMax++;
          }
          if(ID_crit.at(j) == edges.at(i)->Get_name() && p_crit.at(j) <= p_min)
          {
            p_min = p_crit.at(j);
            idxMin = j;
            countMin++;
          }
        }
        if(idxMax >= 0)
        {
          ofsMax << "Pressure limit exceeded in pipe " << edges.at(i)->Get_name() << ", " << countMax << "\t times. Max: ";
          ofsMax << "\t t = "<< t_crit.at(idxMax)<< " s";
          ofsMax << "\t x = "<< x_crit.at(idxMax)<< " m";
          ofsMax << "\t p = "<< 1e-5*p_crit.at(idxMax) << " bar"<< endl;
        }
        else
        {
          //ofsMax << "Pressure limit in pipe " << edges.at(i)->Get_name()  << "\tis NOT exceeded!" << endl;
        }

        if(idxMin >= 0)
        {
          ofsMin << "Pressure limit exceeded in pipe " << edges.at(i)->Get_name() << ", " << countMin << "\t times. Min: ";
          ofsMin << "\t t = "<< t_crit.at(idxMin)<< " s";
          ofsMin << "\t x = "<< x_crit.at(idxMin)<< " m";
          ofsMin << "\t p = "<< 1e-5*p_crit.at(idxMin) << " bar"<< endl;
        }
        else
        {
          //ofsMin << "Pressure limit in pipe " << edges.at(i)->Get_name()  << "\tis NOT exceeded!" << endl;
        }
        
      }
    }
    
    ofsMax.flush();
    ofsMax.close();
    ofsMin.flush();
    ofsMin.close();
    
  }//end of pressure min-max writing

}//end of run function


void PSToolboxRunner::Update_RS(int idx_of_rs, double t_target, bool rs_save)
{
	RS_Solve(idx_of_rs,t_target, rs_save);
	for (int i=0; i< rs_idx_edges.at(idx_of_rs).size(); i++)
		edges.at(rs_idx_edges.at(idx_of_rs).at(i))->UpdateInternal(t_target);

}
