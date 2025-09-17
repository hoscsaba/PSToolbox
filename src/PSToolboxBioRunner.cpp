#include "PSToolboxBioRunner.h"
#include "PSToolboxBaseEdge.h"
#include "Connector.h"
#include "my_tools.h"
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <chrono>

static vector<vector<int> > DEFAULT_VECTOR;

PSToolboxBioRunner::PSToolboxBioRunner(vector<PSToolboxBaseEdge *> &_e,
                                       vector<BioConnector *> &_c, vector<int> &_con_at_edge_start,
                                       vector<int> &_con_at_edge_end, vector<vector<int> > &_rigid_subsystems =
		                                       DEFAULT_VECTOR) {
	edges = _e;
	cons = _c;
	con_at_edge_start = _con_at_edge_start;
	con_at_edge_end = _con_at_edge_end;
	save_data = false;
	folder = "data";
	Node_mul = 1;
	ini_done = false;
	set_limit = false;
	p_limMin = 0.0;
	p_limMax = 0.0;
	write_interval = 1000.0;
	save_interval = 1.0;
	p_start_time = 0.0;
	p_end_time = 1.0e10;
	if (_rigid_subsystems.size() > 0) {
		cout << endl << endl << "ERROR! PSToolboxBioRunner cannot be used with rigid subsystems!!!" << endl;
		cin.get();
	}
	snapshot_fname = "quality";
	snapshot_num = 0;
	save_snapshot = false;
}


void PSToolboxBioRunner::Ini(double dt_target) {
	ini_done = true;
	// Initialization
	for (unsigned int i = 0; i < edges.size(); i++)
		edges.at(i)->Ini(dt_target);
	cout << endl;
	//cin.get();
}

void PSToolboxBioRunner::Ini(double dt_target, VectorXd Cini) {
	for (unsigned int i = 0; i < edges.size(); i++)
		edges.at(i)->Ini(dt_target, Cini);
}


void PSToolboxBioRunner::Run(double t_max) {
	// for elapsed time measurement
	auto start = std::chrono::high_resolution_clock::now();

	cout << endl << endl << "Starting simulation...";


	// Simulation

	double t_global = 0., dt_out = t_max / 100., t_out = -1.e-10;
	double t_next;
	double last_save = -1.e-10;
	double t_save = write_interval;
	int update_idx;
	vector<bool> update_edges(edges.size());
	int STEP = 0;

	double T_TOL = 1.e-10;

	if (save_snapshot) {
		Write_snapshot(t_global);
		snapshot_num++;
	}

	while (t_global < t_max) {
		/*
		 ======================== TIMESTEP PREPARATIONS ===============================
		 */
		STEP++;
		if (DEBUG) {
			cout << endl << "STEP     : " << STEP;
			cout << endl << "t_global : " << t_global;
		}
		if (t_global > t_out) {
			cout << endl << "time: " << convertSecondsToTimeFormat(t_global) << "\t (" << round(t_global / t_max * 100)
					<< "%)";
			t_out += dt_out;
			//cin.get();
		}

		// Find smallest next target time
		fill(update_edges.begin(), update_edges.end(), false);

		update_idx = 0;
		t_next = 1.e15;

		if (DEBUG)
			cout << endl << "Timestep selection ...";
		for (unsigned int i = 0; i < edges.size(); i++) {
			if (!edges.at(i)->is_rigid_element)
				if (edges.at(i)->Get_tnext() < t_next) {
					update_idx = i;
					t_next = edges.at(i)->Get_tnext();

					if (DEBUG)
						printf(
							"\n\t edge #%2d %10s: t=%5.3e, tnext=%5.3e, current tnext=%5.3e",
							i, edges.at(i)->Get_name().c_str(),
							edges.at(i)->Get_t(), edges.at(i)->Get_tnext(),
							t_next);
					//cout << "  <-- t_target ";
				}
		}

		update_edges.at(update_idx) = true;

		// Now find all elements within T_TOL
		for (unsigned int i = 0; i < edges.size(); i++) {
			if (!edges.at(i)->is_rigid_element)
				if (fabs((edges.at(i)->Get_tnext()) - t_next) < T_TOL)
					update_edges.at(i) = true;
		}
		if (DEBUG) {
			cout << endl << endl << "The following edges will be updated: ";
			for (unsigned int i = 0; i < edges.size(); i++) {
				if (update_edges.at(i))
					cout << " " << edges.at(i)->Get_name();
			}
			//cin.get();
		}

		/*
		 ======================== UPDATE EDGES ===============================
		 */

		if (DEBUG)
			cout << endl << "Updating edges... ";

		for (unsigned int i = 0; i < edges.size(); i++) {
			if (update_edges.at(i)) {
				if (DEBUG)
					cout << endl << "\tUpdating edge " << edges.at(i)->name
							<< " ... ->UpdateInternal()";

				edges.at(i)->UpdateInternal(t_next);

				cons.at(con_at_edge_start.at(i))->Update(t_next, i);

				cons.at(con_at_edge_end.at(i))->Update(t_next, i);


				if (DEBUG) {
					cout << endl << " done.";
					//	if (edges.at(i)->name == "801018417") {
					//		cin.get();
					//	}
				}
				if (save_snapshot) {
					if (last_save + save_interval < t_next) {
						Write_snapshot(t_global);
						snapshot_num++;
						last_save += save_interval;
					}
				}
			}
		}
		t_global = t_next;
	}


	// Get the ending time
	auto end = std::chrono::high_resolution_clock::now();

	// Compute total elapsed time in seconds
	auto total_seconds = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();

	// Convert to hours, minutes, and seconds
	int hours = total_seconds / 3600;
	int minutes = (total_seconds % 3600) / 60;
	int seconds = total_seconds % 60;

	// Output the formatted elapsed time
	cout << endl << endl << "Elapsed time: " << hours << " hours, "
			<< minutes << " minutes, "
			<< seconds << " seconds." << std::endl;
} //end of run function

void PSToolboxBioRunner::Save_data_for_ini() {
	cout << endl << endl << "Saving last timestep as ini file to inifile.dat ...";
	FILE *fp;
	fp = fopen("inifile.dat", "w");
	for (auto &edge: edges) {
		fprintf(fp, "%s\t%+5.3e\t%+5.3e\n", edge->name.c_str(), edge->Get_dprop("v_front"),
		        edge->Get_dprop("p_front"));
	}
	fclose(fp);
	cout << " done." << endl;
}


void PSToolboxBioRunner::Write_snapshot(double t_global) {
	string outfname = (snapshot_fname + "_" + std::to_string(snapshot_num) + ".dat");
	cout << endl << " PSToolboxBioRunner::Write_snapshot writing snapshot file " << outfname << " ...";
	FILE *fp;
	fp = fopen(outfname.c_str(), "w");
	fprintf(fp, "Time: %5.3e\n", t_global);
	for (auto &edge: edges) {
		fprintf(fp, "%s\t", edge->name.c_str());
		for (int j = 0; j < edge->Get_iprop("num_of_bio_vars"); j++)
			fprintf(fp, "%+5.3e\t%+5.3e\t", (edge->Get_C_Front())(j), (edge->Get_C_Back())(j));
		fprintf(fp, "\n");
	}
	fclose(fp);
	cout << " done." << endl;
}

static inline std::string ltrim(std::string s) {
	s.erase(s.begin(), std::find_if(s.begin(), s.end(),
	                                [](unsigned char c) { return !std::isspace(c); }));
	return s;
}

void PSToolboxBioRunner::Load_v_conv_from_epanet_export(string fname, double mul) {
	std::cout << "\n\nLoading convective FLOW RATES from " << fname << " ...";

	std::ifstream in(fname);
	if (!in) {
		std::cerr << "\nError: could not open \"" << fname << "\" for reading.\n";
		return;
	}

	// name → edge* map for fast lookups
	std::unordered_map<std::string, PSToolboxBaseEdge *> by_name;
	by_name.reserve(edges.size());
	for (auto *e: edges) {
		by_name.emplace(e->Get_name(), e);
	}

	std::string line;
	std::size_t line_no = 0;
	std::size_t updated = 0, missing = 0, malformed = 0;

	// OPTIONAL: skip the first 3 header lines
	for (int i = 0; i < 3 && std::getline(in, line); ++i) {
		++line_no;
	}

	while (std::getline(in, line)) {
		//cout<<endl<<line;

		++line_no;
		if (line.empty()) continue;

		auto trimmed = ltrim(line);
		if (trimmed.empty()) continue;
		if (trimmed.rfind('#', 0) == 0) continue; // "# comment"
		if (trimmed.rfind("//", 0) == 0) continue; // "// comment"

		std::istringstream iss(trimmed);
		std::string keyword, id;
		double Q_conv = 0.0;

		// Expect:  Pipe <id> <velocity>
		if (!(iss >> keyword >> id >> Q_conv)) {
			std::cerr << "\nWarning: line " << line_no
					<< " malformed; skipping:\n    " << line;
			++malformed;
			continue;
		}
		if (keyword != "Pipe") {
			// Not a pipe line; ignore silently (or count as malformed if you prefer)
			continue;
		}

		auto it = by_name.find(id);
		if (it == by_name.end()) {
			std::cerr << "\nError: edge \"" << id << "\" not found (line " << line_no << ").";
			++missing;
			continue;
		}

		// Q is in l/s
		double set_val = Q_conv / it->second->Get_dprop("A") / 1000.;
		it->second->Set_v_conv(set_val);
		cout << endl << id << "\t: Q=" << Q_conv << " l/s, v=" << set_val << " m/s is set.";
		++updated;
	}

	std::cout << "\nLoaded from " << fname << ". Updated: " << updated
			<< ", Missing: " << missing
			<< ", Malformed: " << malformed << ".\n";
	cout << endl << " Press any key";
	// cin.get();
}

void PSToolboxBioRunner::Load_nodal_demand_from_epanet_export(string fname, double mul) {
	std::cout << "\n\nLoading nodal demands from " << fname << " ...";

	std::ifstream in(fname);
	if (!in) {
		std::cerr << "\nError: could not open \"" << fname << "\" for reading.\n";
		return;
	}

	// name → edge* map for fast lookups
	std::unordered_map<std::string, BioConnector *> by_name;
	by_name.reserve(cons.size());
	for (auto *e: cons) {
		by_name.emplace(e->Get_name(), e);
	}

	std::string line;
	std::size_t line_no = 0;
	std::size_t updated = 0, missing = 0, malformed = 0;

	// OPTIONAL: skip the first 3 header lines
	for (int i = 0; i < 3 && std::getline(in, line); ++i) {
		++line_no;
	}

	while (std::getline(in, line)) {
		cout << endl << line;

		++line_no;
		if (line.empty()) continue;

		auto trimmed = ltrim(line);
		if (trimmed.empty()) continue;
		if (trimmed.rfind('#', 0) == 0) continue; // "# comment"
		if (trimmed.rfind("//", 0) == 0) continue; // "// comment"

		std::istringstream iss(trimmed);
		std::string keyword, id;
		double demand = 0.0;

		// Expect:  Pipe <id> <velocity>
		if (!(iss >> keyword >> id >> demand)) {
			std::cerr << "\nWarning: line " << line_no
					<< " malformed; skipping:\n    " << line;
			++malformed;
			continue;
		}
		if (keyword != "Junc") {
			// Not a pipe line; ignore silently (or count as malformed if you prefer)
			continue;
		}

		auto it = by_name.find(id);
		if (it == by_name.end()) {
			std::cerr << "\nError: node \"" << id << "\" not found (line " << line_no << ").";
			++missing;
			continue;
		}

		// Q is in l/s
		double set_val = demand; // comes in lps, but needed in kg/s
		it->second->Set_demand(set_val);
		cout << endl << id << "\t: demand=" << demand << " l/s is set.";
		++updated;
	}

	std::cout << "\nLoaded from " << fname << ". Updated: " << updated
			<< ", Missing: " << missing
			<< ", Malformed: " << malformed << ".\n";
	cout << endl << " Press any key";
	// cin.get();
}


void PSToolboxBioRunner::Load_v_conv_from_file(string fname, double mul) {
	cout << endl << endl << "Loading convective velocities from " << fname << " ...";

	std::ifstream file(fname); // Open the file
	if (!file) {
		std::cerr << "Error opening file!" << std::endl;
		exit(-1);
	}
	std::ifstream infile(fname);
	if (!infile) {
		std::cerr << "Error: could not open “" << fname << "” for reading.\n";
		return;
	}

	std::string line;
	while (std::getline(infile, line)) {
		if (line.empty()) continue;

		std::istringstream iss(line);
		std::string id_str;
		double v_conv_val, ignore_col3;

		// Read the three columns:
		//   1) id_str (e.g. "p2" or "p4")
		//   2) v_conv_val (e.g. "+4.073e-01")
		//   3) ignore_col3 (e.g. "+8.857e+05")
		if (!(iss >> id_str >> v_conv_val >> ignore_col3)) {
			// malformed line—skip or warn
			std::cerr << "Warning: skipping malformed line:\n    " << line << "\n";
			continue;
		}

		// Find the Edge in edges whose ID matches id_str
		bool is_found = false;

		for (int i = 0; i < edges.size(); i++) {
			if (edges.at(i)->Get_name() == id_str) {
				edges.at(i)->Set_v_conv(limitAbs(v_conv_val, 0.5) * mul);
				is_found = true;
				cout << endl << "\t ID: " << id_str << "\t V_CONV: " << v_conv_val;
				break;
			}
		}

		if (!is_found) {
			std::cerr << "Error: edge " << id_str
					<< " not found;\n";
			cin.get();
		}
	}
	infile.close();
}

void PSToolboxBioRunner::Ini_from_file(string fname, double dt_target) {
	cout << endl << endl << "Loading initial conditions from " << fname << " ...";

	std::ifstream file(fname); // Open the file
	if (!file) {
		std::cerr << "Error opening file!" << std::endl;
		exit(-1);
	}

	std::string line;
	int lineNumber = 0;
	string str;
	double num1, num2;
	while (std::getline(file, line)) {
		lineNumber++;
		std::istringstream iss(line); // Use a string stream to extract numbers
		if (iss >> str >> num1 >> num2) {
			edges.at(lineNumber - 1)->Ini(num1, num2, dt_target);
			//cout<<endl<<"edge "<<edges.at(lineNumber-1)->name<<" v_ini="<<num1<<", p_front="<<num2;
		} else {
			std::cerr << "Error parsing line " << lineNumber << ": " << line << std::endl;
		}
	}
	file.close(); // Close the file
	cout << endl << "Done.";
}

string PSToolboxBioRunner::convertSecondsToTimeFormat(double totalSeconds) {
	int totalSecInt = static_cast<int>(std::floor(totalSeconds)); // truncate fractional part
	int days = totalSecInt / 86400;
	int hours = (totalSecInt % 86400) / 3600;
	int minutes = (totalSecInt % 3600) / 60;

	std::ostringstream oss;
	oss << days << " days "
			<< std::setw(2) << std::setfill('0') << hours << " hours "
			<< std::setw(2) << std::setfill('0') << minutes << " mins ";

	return oss.str();
}

double PSToolboxBioRunner::limitAbs(double value, double maxAbs) {
	if (std::abs(value) > maxAbs) {
		return (value > 0 ? 1 : -1) * maxAbs;
	}
	return value;
}


void PSToolboxBioRunner::build_statistics() {
	double sum_of_inflow = 0.;
	double sum_of_demand = 0.;
	for (auto &con: cons) {
		if (con->demand < 0) {
			sum_of_inflow += con->demand;
			cout << endl << "INFLOW: " << con->name << ", Q=" << con->demand << " l/s";
		} else {
			sum_of_demand += con->demand;
		}
	}
	cout << endl << endl << "sum of inflow: " << -sum_of_inflow << "l/s = " << -sum_of_inflow * 3600 / 1000. << "m3/h";
	cout << endl << "sum of demand: " << sum_of_demand << "l/s = " << sum_of_demand * 3600. / 1000. << "m3/h";

	FILE *outFile = fopen("edges.txt", "w");
	if (!outFile) {
		std::perror("SToolboxBioRunner::build_statistics()\nError opening file"); // prints system error
		cin.get();
	} else
		fprintf(outFile,
		        "name; L (m);  D (m); v_conv (m/s); t_conv (h); t_Monod_2x (h); t_mort_2x (h); t_diff_2x (h) \n");

	for (auto &edge: edges) {
		double v_conv = fabs(edge->Get_dprop("v_conv"));
		double t_conv = edge->Get_dprop("L") / v_conv / 3600.;
		double t_Monod = log(2.) / edge->Get_dprop("bio_mumax") / 3600.;
		double t_mort = log(2.) / edge->Get_dprop("bio_kmort") / 3600.;
		double kfs = edge->Get_dprop("bio_kfs");
		double rh = edge->Get_dprop("D") / 4.;
		double t_diff = 2. / (kfs / rh) / 3600.;
		fprintf(outFile, "%s; %5.3e; %5.3e; %5.3e; %5.3e; %5.3e; %5.3e; %5.3e;\n",
		        edge->Get_name().c_str(), edge->Get_dprop("L"), edge->Get_dprop("D"), v_conv,
		        t_conv, t_Monod, t_mort, t_diff);
		// if (edge->Get_name() == "800014807") {
		// 	cout << endl << "rh =" << edge->Get_dprop("D") / 4.;
		// 	cout << endl << "kfs=" << edge->Get_dprop("bio_kfs");
		// 	cout << endl << "t_diff=" << t_diff;
		// }
	}
	fclose(outFile);
}

