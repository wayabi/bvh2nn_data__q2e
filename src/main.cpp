#include <stdio.h>
#include <math.h>
#include <vector>
#include <deque>
#include <boost/tuple/tuple.hpp>
#include "bvh.h"
#include "qua.h"
#include "thr.h"
#include "rot2.h"
#include "campol.h"
#include "Util.h"
#include "a.h"

typedef boost::math::quaternion<double> Q;

#include <vector>

using namespace std;
using namespace boost;

double get_rand()
{
	return (rand()%1000)/1000.0;
}

int viewer()
{
	BVH bvh;
	//bvh.load("./data/yzx");
	//bvh.load("./data/a_4Char00.bvh");
	bvh.load("./data/scene1.bvh");
	//bvh.load("./data/zai");
	//bvh.load("./data/01_01.bvh");
	ROT2* r = ROT2::make_bone(&bvh);
	double resize0 = r->normalize_height();
	double resize1 = 1.0;
	r->multiply_len(resize1);
	r->update_pos();
	for(int i=0;i<100;++i){
		vector<THR> data = ROT2::get_frame(bvh, i*50);
		r->set_serialized_angle(data);
		r->p_ = r->p_ * resize0*resize1;
		//r->p_ = THR(0, 0, 0);
		double rot = 0;//r->except_y_rotation(r->search("Spine")->q_base_al_cl_.q());
		{
			double x, y, z;
			qua::q2e(r->q_al_cl_.q(), y, x, z, qua::RotSeq::yxz);
			printf("except_y:%f, %f, %f\n", x, y, z);
		}
		//Q qq = qua::e2q(rot, 0, 0, qua::RotSeq::yxz);
		//r->q_al_cl_ = THR(qq*r->q_al_cl_.q());
		//printf("rot_=%f\n", rot-rot_init);

		THR rot_root = r->q_al_cl_;
		r->q_al_cl_ = THR(1, 0, 0, 0);
		ROT2::refer_parent_angle(r);
		r->q_al_cl_ = rot_root;
		r->update_pos();
		//r->print();
		char buf[80];
		sprintf(buf, "./out/%03d.bmp", i);
		//THR dest_camera = r->p_;
		THR dest_camera(0, 0, 0);
		Q q = qua::e2q(0, i*M_PI/50, 0, qua::RotSeq::xyz);
		THR pos_camera = THR(0, 0.4, 2);
		//pos_camera = THR(q*pos_camera.q()/q);
		THR dir_camera = (dest_camera-pos_camera).normalize();
		CamPol::simple_draw(r, buf, pos_camera, dir_camera, NULL);
	}
	delete(r);
	return 0;
}

int t1()
{
	printf("#test Quaternion -> ExponentialMap -> Quaternion\n");
	THR from = THR(get_rand()-0.5, get_rand()-0.5, get_rand()-0.5).normalize();
	THR to = THR(get_rand()-0.5, get_rand()-0.5, get_rand()-0.5).normalize();
	from.print("from");
	to.print("to");
	printf("dot:%f\n", from.dot(to));
	THR q1 = THR(qua::get_quaternion_from_vector(from, to));
	q1.print("q1");
	THR em = qua::q2em(q1.q());
	em.print("em");
	THR q2 = THR(qua::em2q(em));
	q2.print("q2");
	double e = (q1-q2).get_size();
	printf("e=%f\n", e);
	return 0;
}

int t2()
{
	printf("#test interpolation slerp and ExponentialMap\n");
	THR from1 = THR(get_rand()-0.5, get_rand()-0.5, get_rand()-0.5).normalize();
	THR to1 = THR(get_rand()-0.5, get_rand()-0.5, get_rand()-0.5).normalize();
	THR from2 = THR(get_rand()-0.5, get_rand()-0.5, get_rand()-0.5).normalize();
	THR to2 = THR(get_rand()-0.5, get_rand()-0.5, get_rand()-0.5).normalize();
	THR q1 = THR(qua::get_quaternion_from_vector(from1, to1));
	THR q2 = THR(qua::get_quaternion_from_vector(from2, to2));
	THR em1 = qua::q2em(q1.q());
	THR em2 = qua::q2em(q2.q());

	double e = 0.0;
	for(int i=0;i<11;++i){
		double t = i*0.1;
		THR q_interp1 = THR(qua::slerp(q1.q(), q2.q(), t));
		q_interp1.print("q_interp1");
		THR q_interp2 = THR(qua::em2q(em1*(1-t)+em2*t));
		q_interp2.print("q_interp2");
		e += (q_interp1-q_interp2).get_size();
	}
	printf("e=%f\n\n", e);
	return 0;
}

int t3()
{
	double E = 0.0001;
	//printf("#test ExponentialMap manifold");
	THR from = THR(get_rand()-0.5, get_rand()-0.5, get_rand()-0.5).normalize();
	THR to = THR(get_rand()-0.5, get_rand()-0.5, get_rand()-0.5).normalize();
	THR q1 = THR(qua::get_quaternion_from_vector(from, to));
	THR em = qua::q2em(q1.q());
	em.print("em");
	THR q2 = THR(qua::em2q(em));
	double e = (q1-q2).get_size();
	if(e*e > E) printf("error e=%f\n", e);
	return 0;
}

//decompose root dir
int t4()
{
	printf("#test decompose\n");
	THR from = THR(get_rand()-0.5, get_rand()-0.5, get_rand()-0.5).normalize();
	THR to = THR(get_rand()-0.5, get_rand()-0.5, get_rand()-0.5).normalize();
	from.print("from");
	to.print("to");
	THR q1 = THR(qua::get_quaternion_from_vector(from, to));

	THR origin(1, 0, 0);
	THR dir_swing = THR(q1.q()*origin.q()/q1.q());
	THR swing = THR(qua::get_quaternion_from_vector(dir_swing, origin));
	double v1 = 0;
	double v2 = 0;
	double v3 = 0;
	qua::q2e(swing.q(), v1, v2, v3, qua::RotSeq::xzy);
	THR q2 = THR(qua::e2q(0, 0, v3, qua::RotSeq::xzy));
	
	THR q3 = q1.q()/q2.q();
	THR q4 = q3.q()*q2.q();
	q2.print("decompose(y)");
	q3.print("decompose(x, z)");
	q1.print("origin ");
	q4.print("de->com");

	THR dir_q1 = THR(q1.q()*origin.q()/q1.q());
	THR dir_q2 = THR(q2.q()*origin.q()/q2.q());
	printf("tan(q1_y_axis) = %f\n", tan(dir_q1.x_/dir_q1.z_));
	printf("tan(q2_y_axis) = %f\n", tan(dir_q2.x_/dir_q2.z_));
	return 0;
}

//decompose y_rotation
int t5()
{
	BVH bvh;
	bvh.load("./data/01_01.bvh");
	//bvh.load("./data/yzx");
	//bvh.load("./data/zai");
	ROT2* r = ROT2::make_bone(&bvh);
	r->update_pos();
	double resize0 = r->normalize_height();
	double resize1 = 50;
	r->multiply_len(resize1);
	for(int i=0;i<100;++i){
		vector<THR> data = ROT2::get_frame(bvh, i*10);
		r->set_serialized_angle(data);
		auto except_y = r->except_y_rotation();
		double y_rot = except_y.first;
		r->q_al_cl_ = except_y.second;
		r->multiply_pos(resize0*resize1);
		printf("y_rot = %f\n", y_rot);
		r->update_pos();
		//r->print();
		char buf[80];
		sprintf(buf, "./out/%03d.bmp", i);
		//THR dest_camera = r->p_;
		THR dest_camera(0, 0, 0);
		Q q = qua::e2q(0, 0*i*M_PI/50, 0, qua::RotSeq::xyz);
		THR pos_camera = THR(0, 35, 100);
		pos_camera = THR(q*pos_camera.q()/q);
		THR dir_camera = (dest_camera-pos_camera).normalize();
		CamPol::simple_draw(r, buf, pos_camera, dir_camera, NULL);
	}
	delete(r);
	return 0;
}

//rots load
int t6()
{
	BVH bvh;
	//bvh.load("./data/01_01.bvh");
	//bvh.load("./data/yzx");
	//bvh.load("./data/zai");
	//bvh.load("./data/a_4Char00.bvh");
	bvh.load("./data/scene1.bvh");

	std::vector<ROT2*> rots;
	std::vector<THR> pos_diff;
	std::vector<double> y_rot_diff;
	std::vector<double> y_root;
	THR init_pos(0, 0, 0);
	THR last_pos(0, 0, 0);
	double init_y_rot = 0;
	double last_y_rot = 0;

	int i=0;
	for(auto ite = bvh.motion_.begin();ite != bvh.motion_.end();++ite, ++i){
		ROT2* r = ROT2::make_bone(&bvh);
		double resize0 = r->normalize_height();
		r->update_pos();
		vector<THR> data = ROT2::get_frame(bvh, i);
		r->set_serialized_angle(data);
		auto except_y = r->except_y_rotation();
		double y_rot = except_y.first;
		r->q_al_cl_ = except_y.second;
		if(i==0){
			last_pos = r->p_;
			init_pos = r->p_;
			last_y_rot = y_rot;
			init_y_rot = y_rot;
		}
		//r->p_ = (r->p_-init_pos)*resize0;
		//r->p_ = THR(q*r->p_.q()/q);
		//r->p_ = THR(0, 0, 0);
		//Q q = qua::get_quaternion_from_axis(0, 1, 0, -y_rotate);
		Q q = qua::e2q(-init_y_rot, 0, 0, qua::RotSeq::yxz);
		//THR(q).print("q");
		//THR p2 = THR(q*r->p_.q()/q);
		//THR init_pos2 = THR(q*init_pos.q()/q);
		//r->p_ = (p2-init_pos2)*resize0;
		THR pos = (r->p_-last_pos)*resize0;
		//THR pos = (r->p_-init_pos)*resize0;
		pos = THR(q*pos.q()/q);
		pos_diff.push_back(pos);
		//THR temp = (r->p_-temp_init_pos)*resize0;
		//temp = THR(q*temp.q()/q);
		last_pos = r->p_;
		r->p_ = THR(0, 0, 0);
		r->update_pos();
		vector<double> minmax = r->get_minmax_pos();
		y_root.push_back(-minmax.at(2));
		printf("min=%f, max=%f\n", minmax.at(2), minmax.at(3));
		//r->p_ = temp;
		
		y_rot_diff.push_back(y_rot-last_y_rot);
		last_y_rot = y_rot;

		rots.push_back(r);
	}

	THR pos_sum(0, 0, 0);
	double rot_sum = 0;//init_y_rot;
	int skip = 20;
	for(int i=0;i<100;++i){
		for(int j=0;j<skip;++j){
			if(i == 0) break;
			pos_sum = pos_sum + pos_diff.at((i-1)*skip+j);
			rot_sum = rot_sum + y_rot_diff.at((i-1)*skip+j);
		}
		ROT2* r = rots.at(i*skip);
		r->p_ = r->p_+pos_sum;
		Q q = qua::e2q(rot_sum, 0, 0, qua::RotSeq::yxz);
		r->q_al_cl_ = THR(q*r->q_al_cl_.q());

		THR rot_root = r->q_al_cl_;
		r->q_al_cl_ = THR(1, 0, 0, 0);
		ROT2::refer_parent_angle(r);
		r->q_al_cl_ = rot_root;

		r->update_pos();
		char buf[80];
		sprintf(buf, "./out/%03d.bmp", i);
		//THR dest_camera = r->p_;
		THR dest_camera(0, 0, 0);
		THR pos_camera = THR(0, 1, 3);
		//pos_camera = THR(q*pos_camera.q()/q);
		THR dir_camera = (dest_camera-pos_camera).normalize();
		CamPol::simple_draw(r, buf, pos_camera, dir_camera, NULL);
	}

	for(auto ite = rots.begin();ite != rots.end();++ite){
		delete(*ite);
	}
	return 0;
}

//decompose
int t7()
{
	THR from1 = THR(get_rand()-0.5, get_rand()-0.5, get_rand()-0.5).normalize();
	THR to1 = THR(get_rand()-0.5, get_rand()-0.5, get_rand()-0.5).normalize();
	THR from2 = THR(get_rand()-0.5, get_rand()-0.5, get_rand()-0.5).normalize();
	THR to2 = THR(get_rand()-0.5, get_rand()-0.5, get_rand()-0.5).normalize();
	THR q11 = THR(qua::get_quaternion_from_vector(to1, from1));
	THR q12 = THR(qua::get_quaternion_from_vector(to2, from2));
	THR q = q11.q()*q12.q();

	double v1 = 0;
	double v2 = 0;
	double v3 = 0;
	qua::q2e(q.q(), v1, v2, v3, qua::RotSeq::yxz);
	THR q2 = THR(qua::e2q(v1, 0, 0, qua::RotSeq::yxz));
	THR q3 = THR(qua::e2q(0, v2, v3, qua::RotSeq::yxz));

	//printf("1nearly_equal() = %d\n", q.nearly_equal(q3.q()*q2.q()));
	bool flag = q.nearly_equal(q2.q()*q3.q());
	printf("2nearly_equal() = %d\n", flag);
	if(!flag){
		q.print("q1");
		THR(q2.q()*q3.q()).print("q2");
	}
	
	return 0;
}

std::pair<THR, THR> get_stick(ROT2* r, double len)
{
	THR origin(len, 0, 0);
	THR q_root = r->q_al_cl_;
	THR q_stick = THR(qua::e2q(-M_PI/2, M_PI/2, 0, qua::RotSeq::yzx));
	//THR q_al_cw = THR(q_root.q()*q_stick.q()/q_root.q());
	THR q_al_cw = q_stick;
	THR q = THR(q_root.q()*q_al_cw.q());
	//THR q = q_root;
	THR dir = q.q()*origin.q()/q.q();

	return std::make_pair(r->p_-dir, q);
}

//auto stick add
int t8()
{
	double len_stick = 40;

	BVH bvh;
	//bvh.load("./data/01_01.bvh");
	//bvh.load("./data/yzx");
	//bvh.load("./data/zai");
	bvh.load("./data/scene1.bvh");
	//bvh.load("./data/a_4Char00.bvh");

	ROT2* r = ROT2::make_bone(&bvh);
	double resize0 = r->normalize_height();
	double resize1 = 80;
	r->multiply_len(resize1);
	r->update_pos();
	for(int i=0;i<100;++i){
		vector<THR> data = ROT2::get_frame(bvh, i*100);
		r->set_serialized_angle(data);
		r->p_ = r->p_*resize0*resize1;
		r->update_pos();
		//r->print();
		char buf[80];
		sprintf(buf, "./out/%03d.bmp", i);
		//THR dest_camera = r->p_;
		THR dest_camera(0, 100, 0);
		dest_camera = r->p_;
		Q q = qua::e2q(0, i*M_PI/50, 0, qua::RotSeq::xyz);
		THR pos_camera = THR(160, 120, 0);
		//pos_camera = THR(q*pos_camera.q()/q);
		THR dir_camera = (dest_camera-pos_camera).normalize();

		auto stick = get_stick(r, len_stick);
		THR pos_stick = stick.first;
		THR q_stick = stick.second;
		vector<vector<PolColor> > bones;
		CamPol cp;
		bones.push_back(cp.make_bone(pos_stick, q_stick, len_stick));
		CamPol::simple_draw(r, buf, pos_camera, dir_camera, &bones);
	}
	delete(r);
	return 0;
}

void write_nn_input(FILE* f, std::vector<THR>& em_input, double y_root, THR& pos_stick_f, THR& em_stick_f)
{
	for(auto ite = em_input.begin();ite != em_input.end();++ite){
		fprintf(f, "%f,%f,%f,", ite->x_, ite->y_, ite->z_);
	}
	fprintf(f, "%f,%f,%f,%f,%f,%f,%f\n", y_root, pos_stick_f.x_, pos_stick_f.y_, pos_stick_f.z_, em_stick_f.x_, em_stick_f.y_, em_stick_f.z_);
}

void write_nn_input_v(FILE* f, std::vector<THR>& em_input, double y_root, THR& pos_stick_f, THR& em_stick_f, vector<THR>& v_part)
{
	for(auto ite = em_input.begin();ite != em_input.end();++ite){
		fprintf(f, "%f,%f,%f,", ite->x_, ite->y_, ite->z_);
	}
	fprintf(f, "%f,%f,%f,%f,%f,%f,%f", y_root, pos_stick_f.x_, pos_stick_f.y_, pos_stick_f.z_, em_stick_f.x_, em_stick_f.y_, em_stick_f.z_);

	int size = 47;
	int index_unity_bone_exist[] = {
		/*0,*/7,9,11,36,37,38,39,40,41,42,44,45,46,48,49,50,52,53,54,56,57,58,13,14,15,16,17,18,19,21,22,23,25,26,27,29,30,31,33,34,35,4,5,6,1,2,3
	};
	for(int i=0;i<size;++i){
		THR t = v_part.at(index_unity_bone_exist[i]);
		fprintf(f, ",%f,%f,%f", t.x_, t.y_, t.z_);
	}
	fprintf(f, "\n");
}

void write_nn_output(FILE* f, std::vector<THR>& em_output, THR& pos_output, double y_rot_output)
{
	for(auto ite = em_output.begin();ite != em_output.end();++ite){
		fprintf(f, "%f,%f,%f,", ite->x_, ite->y_, ite->z_);
	}
	fprintf(f, "%f,%f,%f,%f\n", pos_output.x_, pos_output.y_, pos_output.z_, y_rot_output);
}

int make_nn_input()
{
	BVH bvh;
	//bvh.load("./data/01_01.bvh");
	//bvh.load("./data/yzx");
	//bvh.load("./data/zai");
	//bvh.load("./data/a_4Char00.bvh");
	//bvh.load("./data/a");
	bvh.load("./data/scene1.bvh");

	FILE* f_nn_input = NULL;
	FILE* f_nn_output = NULL;
	if((f_nn_input = fopen("f_nn_input.txt", "w")) == NULL){
		printf("file open error1.\n");
		return 1;
	}
	if((f_nn_output = fopen("f_nn_output.txt", "w")) == NULL){
		printf("file open error2.\n");
		fclose(f_nn_input);
		return 1;
	}

	vector<int> frame_feature_stick;
	//frame_feature_stick.push_back(0);
	//frame_feature_stick.push_back(1);
/*
	frame_feature_stick.push_back(60);
	frame_feature_stick.push_back(120);
	frame_feature_stick.push_back(180);
	frame_feature_stick.push_back(240);
*/
	//for(int i=0;i<10;++i){
	//	frame_feature_stick.push_back((i+1)*30);
	//}
	frame_feature_stick.push_back(1);
	int frame_feature_max = frame_feature_stick.at(frame_feature_stick.size()-1);
	double len_stick = 0.3;

	int frame_skip = 0;
	int frame_next = 1;
	int frame_start = 0;
	int frame_end = bvh.motion_.size()/(frame_skip+1);
	//int frame_end = 241;

	vector<ROT2*> rots;
	vector<THR> pos_diff;
	vector<double> y_rot_diff;
	vector<double> y_root;
	vector<pair<THR, THR> > sticks;
	THR init_pos(0, 0, 0);
	THR last_pos(0, 0, 0);
	double init_y_rot = 0;
	double last_y_rot = 0;

	int count_skip = 0;
	int i=0;
	for(auto ite = bvh.motion_.begin();ite != bvh.motion_.end();++ite, ++i, ++count_skip){
		if(count_skip >= frame_skip){
			count_skip = 0;
		}else{
			continue;
		}
		ROT2* r = ROT2::make_bone(&bvh);
		double resize0 = r->normalize_height();
		r->update_pos();
		vector<THR> data = ROT2::get_frame(bvh, i);
		r->set_serialized_angle(data);
		auto except_y = r->except_y_rotation(last_y_rot);
		double y_rot = except_y.first;
		r->q_al_cl_ = except_y.second;
{
		double x = 0;
		double y = 0;
		double z = 0;
		qua::q2e(r->q_al_cl_.q(), y, x, z, qua::RotSeq::yxz);
		printf("y_extract(%f, %f, %f)\n", x, y, z);
}

		if(rots.size() == 0){
			last_pos = r->p_;
			init_pos = r->p_;
			last_y_rot = y_rot;
			init_y_rot = y_rot;
		}
		Q q = qua::e2q(-y_rot, 0, 0, qua::RotSeq::yxz);
		THR pos = (r->p_-last_pos)*resize0;
		pos = THR(q*pos.q()/q);
		pos_diff.push_back(pos);
		last_pos = r->p_;
		r->p_ = THR(0, 0, 0);
		r->update_pos();
		sticks.push_back(get_stick(r, len_stick));
		vector<double> minmax = r->get_minmax_pos();
		y_root.push_back(-minmax.at(2));
		//printf("min=%f, max=%f\n", minmax.at(2), minmax.at(3));
		
		y_rot_diff.push_back(y_rot-last_y_rot);
		last_y_rot = y_rot;

		rots.push_back(r);
	}

	//###nn
	for(int i=frame_start;i<frame_end-frame_feature_max;++i){
		for(int j=0;j<(int)frame_feature_stick.size();++j){
			int frame_feature = frame_feature_stick.at(j);
			THR pos_sum(0, 0, 0);
			double rot_sum = 0;
			for(int k=i;k<=i+frame_feature;++k){
				pos_sum = pos_sum+pos_diff.at(k);
				rot_sum = rot_sum+y_rot_diff.at(k);
			}
			Q q = qua::e2q(rot_sum, 0, 0, qua::RotSeq::yxz);
			THR pos_stick_f = THR(q*sticks.at(i+frame_feature).first.q()/q)+pos_sum;
			THR em_stick_f = qua::q2em(q*sticks.at(i+frame_feature).second.q());
			
			vector<THR> angle_input = rots.at(i)->get_serialized_angle();
			//vector<THR> angle_input = rots.at(i)->get_serialized_angle_al_cw();
			vector<THR> em_input;
			for(auto ite = angle_input.begin();ite != angle_input.end();++ite){
				if(ite != angle_input.begin()){
					//skip first pos_data
					em_input.push_back(qua::q2em(ite->q()));
					double y = 0;
					double x = 0;
					double z = 0;
					qua::q2e(ite->q(), y, x, z, qua::RotSeq::yxz);
					//printf("y_extract3(%f, %f, %f)\n", x*180/M_PI, y*180/M_PI, z*180/M_PI);
				}
			}
			int ii=0;
			//for unity humanoid. Needs al_cw.
			//vector<THR> angle_output = rots.at(i+frame_next)->get_serialized_angle_al_cw();
			vector<THR> angle_output = rots.at(i+frame_next)->get_serialized_angle();
			vector<THR> em_output;
			for(auto ite = angle_output.begin();ite != angle_output.end();++ite, ++ii){
				if(ite != angle_output.begin()){
					double y = 0;
					double x = 0;
					double z = 0;
					qua::q2e(ite->q(), y, x, z, qua::RotSeq::yxz);
					if(ii==1){
						printf("y_extract2(%f, %f, %f)\n", x*180/M_PI, y*180/M_PI, z*180/M_PI);
						ite->print();
					}
						
					//*ite = THR(qua::e2q::e2q(y, x, z, qua::RotSeq::zxy));
					//skip first pos_data
					THR t = *ite;
					//t.x_ *= -1;
					//t.w_ *= -1;
					THR em = qua::q2em(t.q());
					THR qq = qua::em2q(em);
					double error = (t-qq).get_size();
					if(error > 0.1){
						printf("error!\n");
						t.print("1");
						qq.print("2");
					}
					//em.y_ *= -1;
					//em.z_ *= -1;
					em_output.push_back(em);
				}
			}
			double y_root_input = y_root.at(i);
			THR pos_output = THR(0, 0, 0);
			double y_rot_output = 0;
			for(int k=i;k<=i+frame_next;++k){
				pos_output = pos_output+pos_diff.at(k);
				y_rot_output = y_rot_output+y_rot_diff.at(k);
			}
			y_rot_output = qua::convert_single_pi(y_rot_output);

			write_nn_input(f_nn_input, em_input, y_root_input, pos_stick_f, em_stick_f);
			printf("aaa %f\n", y_rot_output);
			write_nn_output(f_nn_output, em_output, pos_output, y_rot_output);
			printf("%d,%d,%d\n", i, j, frame_feature);
		}
	}

	//draw
	if(false){
	THR pos_sum(0, 0, 0);
	double rot_sum = 0;//init_y_rot;
	for(int i=0;i<100;++i){
		for(int j=0;j<=frame_skip;++j){
			if(i == 0) break;
			pos_sum = pos_sum + pos_diff.at((i-1)*(frame_skip+1)+j);
			rot_sum = rot_sum + y_rot_diff.at((i-1)*(frame_skip+1)+j);
		}
		ROT2* r = rots.at(i);
		r->p_ = r->p_+pos_sum;
		Q q = qua::e2q(rot_sum, 0, 0, qua::RotSeq::yxz);
		//r->q_al_cl_ = THR(q*r->q_al_cl_.q());
		r->update_pos();
		char buf[80];
		sprintf(buf, "./out/%03d.bmp", i);
		THR dest_camera(0, 0, 0);
		THR pos_camera = THR(3, 1, 0);
		THR dir_camera = (dest_camera-pos_camera).normalize();

		//stick
		vector<vector<PolColor> > stick_bone;
		CamPol cp;
		THR pos_stick = THR(q*sticks.at(i*10).first.q()/q)+r->p_;
		THR q_stick = THR(q*sticks.at(i*10).second.q());//THR(q*sticks.at(i*10).second.q());
		stick_bone.push_back(cp.make_bone(pos_stick, q_stick, len_stick));

		CamPol::simple_draw(r, buf, pos_camera, dir_camera, &stick_bone);
	}
	}

	for(auto ite = rots.begin();ite != rots.end();++ite){
		delete(*ite);
	}

	fclose(f_nn_input);
	fclose(f_nn_output);
	return 0;
}

void t_rotation_order()
{
	THR origin(1, 0, 0);

	Q q = qua::e2q(0, M_PI/2, M_PI/4, qua::RotSeq::xyz);
	THR(q*origin.q()/q).print();
}

void tmp()
{
	THR a(0.2, -90, -90);
	THR b(80.9, 0, 317.1);
	a = a*M_PI/180;
	b = b*M_PI/180;
	Q qa = qua::e2q(a.y_, a.x_, a.z_, qua::RotSeq::yxz);
	Q qb = qua::e2q(b.y_, b.x_, b.z_, qua::RotSeq::yxz);
	double x, y, z;
	qua::q2e(qb*qa, y, x, z, qua::RotSeq::yxz);
	THR(qa).print("qa");
	THR(qb).print("qb");
	printf("qb*qa:%f, %f, %f\n", 180*x/M_PI, 180*y/M_PI, 180*z/M_PI);
	
}

vector<vector<THR> > load_em2data()
{
	vector<vector<THR> > ret;
	string path_file = "f_nn_output.txt";
	auto sss = Util::load_csv(path_file.c_str());
	for(auto ite = sss.begin();ite != sss.end();++ite){
		vector<THR> qq;
		qq.push_back(THR(0, 0, 0));
		for(int i=0;i<(ite->size()-4)/3;++i){
			float x = (float)atof(ite->at(i*3+0).c_str());
			float y = (float)atof(ite->at(i*3+1).c_str());
			float z = (float)atof(ite->at(i*3+2).c_str());
			THR t(x, y, z);
			THR q = qua::em2q(t);
			qq.push_back(q);
		}
		float x = atof(ite->at(ite->size()-4).c_str());
		float y = atof(ite->at(ite->size()-3).c_str());
		float z = atof(ite->at(ite->size()-2).c_str());
		float y_rot = atof(ite->at(ite->size()-1).c_str());
		qq.push_back(THR(y_rot, x, y, z));
		ret.push_back(qq);
	}
	return ret;
}
//load em nn_output
int t9()
{
	BVH bvh;
	//bvh.load("./data/yzx");
	//bvh.load("./data/a_4Char00.bvh");
	bvh.load("./data/scene1.bvh");
	//bvh.load("./data/zai");
	//bvh.load("./data/01_01.bvh");
	ROT2* r = ROT2::make_bone(&bvh);
	double resize0 = r->normalize_height();
	double resize1 = 1.0;
	r->multiply_len(resize1);
	r->update_pos();
	double rot_init = 0;
	auto data_all = load_em2data();
	printf("size = %d\n", data_all.size());
	THR pos(0, 0, 0);
	for(int i=0;i<100;++i){
		vector<THR> data = data_all.at(i*1);
		r->set_serialized_angle(data);
		double y_rot = data.at(60).w_;
		Q q_y_rot = qua::e2q(-y_rot, 0, 0, qua::RotSeq::yxz);
		pos = pos+THR(q_y_rot*THR(data.at(60)).q()/q_y_rot)*resize1;
		r->p_ = pos;
		{
			double x, y, z;
			qua::q2e(r->q_al_cl_.q(), y, x, z, qua::RotSeq::yxz);
			printf("except_y:%f, %f, %f\n", x, y, z);
		}
		//Q qq = qua::e2q(rot, 0, 0, qua::RotSeq::yxz);
		//r->q_al_cl_ = THR(qq*r->q_al_cl_.q());
		//printf("rot_=%f\n", rot-rot_init);

		THR rot_root = r->q_al_cl_;
		r->q_al_cl_ = THR(1, 0, 0, 0);
		ROT2::refer_parent_angle(r);
		r->q_al_cl_ = rot_root;
		r->update_pos();
		if(i== 29 || i == 30 || i==31){
			printf("#### %d\n", i);
			r->print();
		}
		//r->print();
		char buf[80];
		sprintf(buf, "./out/%03d.bmp", i);
		//THR dest_camera = r->p_;
		THR dest_camera(0, 0, 0);
		Q q = qua::e2q(0, i*M_PI/50, 0, qua::RotSeq::xyz);
		THR pos_camera = THR(0, 0.4, 2);
		//pos_camera = THR(q*pos_camera.q()/q);
		THR dir_camera = (dest_camera-pos_camera).normalize();
		CamPol::simple_draw(r, buf, pos_camera, dir_camera, NULL);
	}
	delete(r);
	return 0;
}

void t10()
{
	THR a(180, 60, 180);
	THR b(-180, 60, 180);
	a = a*M_PI/180;
	b = b*M_PI/180;
	a.print("a");
	b.print("b");
	Q qa = qua::e2q(a.y_, a.x_, a.z_, qua::RotSeq::yxz);
	Q qb = qua::e2q(b.y_, b.x_, b.z_, qua::RotSeq::yxz);
	THR(qa).print("qa");
	THR(qb).print("qb");
	THR aa;
	THR bb;
	qua::q2e(qa, aa.x_, aa.y_, aa.z_, qua::RotSeq::yxz);
	qua::q2e(qb, bb.x_, bb.y_, bb.z_, qua::RotSeq::yxz);
	aa.print("euler_qa");
	bb.print("euler_qb");

	printf("\n");

	THR ema = qua::q2em(qa);
	THR emb = qua::q2em(qb);
	ema.print("ema");
	emb.print("emb");
	Q q_ema = qua::em2q(ema);
	Q q_emb = qua::em2q(emb);
	THR(q_ema).print("q_ema");
	THR(q_emb).print("q_emb");
	THR emaa;
	THR embb;
	qua::q2e(q_ema, emaa.x_, emaa.y_, emaa.z_, qua::RotSeq::yxz);
	qua::q2e(q_emb, embb.x_, embb.y_, embb.z_, qua::RotSeq::yxz);
	emaa.print("euler_q_ema");
	embb.print("euler_q_emb");
	printf("qb[x] = %.32f\n", qb.R_component_2());
	printf("q_emb[x] = %.32f\n", q_emb.R_component_2());
}

//except_y_rot
int t11()
{
	BVH bvh;
	//bvh.load("./data/01_01.bvh");
	//bvh.load("./data/yzx");
	//bvh.load("./data/zai");
	bvh.load("./data/a_4Char00.bvh");

	std::vector<ROT2*> rots;
	std::vector<THR> pos_diff;
	std::vector<double> y_rot_diff;
	std::vector<double> y_root;
	THR init_pos(0, 0, 0);
	THR last_pos(0, 0, 0);
	double init_y_rot = 0;
	double last_y_rot = 0;

	vector<THR> qq1;
	vector<THR> pos1;

	int i=0;
	for(auto ite = bvh.motion_.begin();ite != bvh.motion_.end();++ite, ++i){
		ROT2* r = ROT2::make_bone(&bvh);
		double resize0 = r->normalize_height();
		r->update_pos();
		vector<THR> data = ROT2::get_frame(bvh, i);
		r->set_serialized_angle(data);
		qq1.push_back(r->q_al_cl_);
		auto except_y = r->except_y_rotation();
		double y_rot = except_y.first;
		r->q_al_cl_ = except_y.second;
		
		if(i==0){
			last_pos = r->p_;
			init_pos = r->p_;
			last_y_rot = 0;//y_rot;
			init_y_rot = y_rot;
		}
		Q q = qua::e2q(-init_y_rot, 0, 0, qua::RotSeq::yxz);
		THR pos = (r->p_-last_pos)*resize0;
		pos1.push_back(pos);
		pos = THR(q*pos.q()/q);
		pos_diff.push_back(pos);
		last_pos = r->p_;
		r->p_ = THR(0, 0, 0);
		r->update_pos();
		vector<double> minmax = r->get_minmax_pos();
		y_root.push_back(-minmax.at(2));
		//printf("min=%f, max=%f\n", minmax.at(2), minmax.at(3));
		//r->p_ = temp;
		
		y_rot_diff.push_back(y_rot-last_y_rot);
		last_y_rot = y_rot;

		rots.push_back(r);
	}

	//main logic
	int ii = 0;
	double y_rot_sum = 0;
	for(auto ite = rots.begin();ite != rots.end();++ite, ++ii){
		THR q = (*ite)->q_al_cl_;
		double y_rot = y_rot_diff.at(ii);
		y_rot_sum += y_rot;
		Q q_y_rot = qua::e2q(y_rot_sum, 0, 0, qua::RotSeq::yxz);
		THR tq = THR(q_y_rot*q.q());
		THR qq = qq1.at(ii);
		qq.print("qq");
		tq.print("tq");
		printf("q_error = %f\n", (qq-tq).get_size());
		

		THR pos = pos1.at(ii);

		
	}

	for(auto ite = rots.begin();ite != rots.end();++ite){
		delete(*ite);
	}
	return 0;
}

void t12()
{
	THR from = THR(get_rand()-0.5, get_rand()-0.5, get_rand()-0.5).normalize();
	THR to = THR(get_rand()-0.5, get_rand()-0.5, get_rand()-0.5).normalize();
	//from.print("from");
	//to.print("to");
	//printf("dot:%f\n", from.dot(to));
	THR q1 = THR(qua::get_quaternion_from_vector(from, to));

	double y_rot = get_rand()*M_PI*2-M_PI;
	THR q2 = THR(qua::e2q(-y_rot, 0, 0, qua::RotSeq::yxz));
	THR q3 = THR(q2.q()*q1.q());
	THR q4 = THR(qua::e2q(y_rot, 0, 0, qua::RotSeq::yxz)*q3.q());
	//q1.print("q1");
	//q4.print("q4");
	printf("error = %f\n", (q1-q4).get_size());
}

vector<double> get_nn_input(ROT2* r, THR pos_stick, Q q_stick, double y_rot_sum)
{
	Q q = qua::e2q(0, 0, -y_rot_sum, qua::RotSeq::zxy);
	vector<double> data_input;
	vector<THR> data = r->get_serialized_angle();
	for(auto ite = data.begin();ite != data.end();++ite){
		if(ite != data.begin()){
			THR em = qua::q2em(ite->q());
			data_input.push_back(em.x_);
			data_input.push_back(em.y_);
			data_input.push_back(em.z_);
		}
	}
	vector<double> minmax = r->get_minmax_pos();
	double y_root = r->p_.y_-minmax.at(2);
	data_input.push_back(y_root);

	pos_stick = pos_stick - r->p_;
	pos_stick = THR(q*pos_stick.q()/q);
	q_stick = q_stick*q;
	data_input.push_back(pos_stick.x_);
	data_input.push_back(pos_stick.y_);
	data_input.push_back(pos_stick.z_);
	THR em = qua::q2em(q_stick);
	data_input.push_back(em.x_);
	data_input.push_back(em.y_);
	data_input.push_back(em.z_);
	return data_input;
	
}

vector<double> do_nn(vector<double>& data_input)
{
	vector<double> ret;
	string file_input = "./tmp/nn_input";
	string file_output = "./tmp/nn_output";
	FILE* f;
	if((f = fopen(file_input.c_str(), "w")) == NULL){
		return ret;
	}
	for(int i=0;i<data_input.size();++i){
		if(i==0){
			fprintf(f, "%f", data_input.at(i));
		}else{
			fprintf(f, ",%f", data_input.at(i));
		}
	}
	fprintf(f, "\n");
	fclose(f);
	system((string("python /home/ambi/myc/keras_bvh/bvh_nn.py ") + file_input + string(" > ") + file_output).c_str());

	auto csv = Util::load_csv(file_output.c_str());
	for(auto ite = csv.at(0).begin();ite != csv.at(0).end();++ite){
		ret.push_back(atof(ite->c_str()));
	}
	return ret;
}

int t_nn1()
{
	BVH bvh;
	bvh.load("./data/scene1.bvh");
	ROT2* r = ROT2::make_bone(&bvh);
	double resize0 = r->normalize_height();
	double resize1 = 1.0;
	r->multiply_len(resize1);
	vector<THR> data = ROT2::get_frame(bvh, 100);
	r->set_serialized_angle(data);
	r->p_ = THR(0, 0, 0);
	THR pos_stick(0, 0.3, -0.1);
	Q q_stick = qua::e2q(0, 0, M_PI/2, qua::RotSeq::xyz);
	double y_rot_sum = 0;
	auto except_y = r->except_y_rotation();
	double y_rot = except_y.first;
	r->q_al_cl_ = except_y.second;
	y_rot_sum += y_rot;
	r->update_pos();

	for(int i=0;i<100;++i){
		vector<double> data = get_nn_input(r, pos_stick, q_stick, y_rot_sum);
		vector<double> output = do_nn(data);

		vector<THR> rots;
		THR pos = THR(output.at(59*3+0), output.at(59*3+1), output.at(59*3+2));
		y_rot_sum += output.at(59*3+3);
		Q q_y_rot = qua::e2q(0, 0, -y_rot_sum, qua::RotSeq::zxy);
		pos = THR(q_y_rot*pos.q()/q_y_rot);
		rots.push_back(pos);
		for(int j=0;j<59;++j){
			THR em(output.at(j*3+0), output.at(j*3+1), output.at(j*3+2));
			THR q = THR(qua::em2q(em));
			rots.push_back(q);
		}
		r->set_serialized_angle(rots);

		THR rot_old = r->q_al_cw_;
		r->q_al_cl_ = THR(conj(q_y_rot)*r->q_al_cl_.q()*q_y_rot);
		r->update_pos();
		r->q_al_cl_ = rot_old;
		char buf[80];
		sprintf(buf, "./out/%03d.bmp", i);
		THR dest_camera(0, 0, 0);
		Q q = qua::e2q(0, i*M_PI/50, 0, qua::RotSeq::xyz);
		THR pos_camera = THR(0, 0.4, 2);
		THR dir_camera = (dest_camera-pos_camera).normalize();
		CamPol::simple_draw(r, buf, pos_camera, dir_camera, NULL);
	}
	delete(r);
	return 0;
}

int t13()
{
	BVH bvh;
	bvh.load("./data/scene1.bvh");
	ROT2* r = ROT2::make_bone(&bvh);
	double resize0 = r->normalize_height();
	double resize1 = 1.0;
	r->multiply_len(resize1);

	int frame_start = 3000;
	int frame_end = 3300;

	FILE* f_input;
	FILE* f_output;
	f_input = fopen("./out/nn1_input", "w");
	f_output = fopen("./out/nn1_output", "w");
	
	THR pos_last(0, 0, 0);
	double rot_y_last = 0;
	for(int i=frame_start;i<=frame_end;i+=5){
		double weight = (i-frame_start)/(double)(frame_end-frame_start);
		vector<THR> data = ROT2::get_frame(bvh, i);
		r->set_serialized_angle(data);
		auto except_y = r->except_y_rotation();
		double rot_y = except_y.first;
		r->q_al_cl_ = except_y.second;
		THR pos_diff;
		double rot_y_diff = 0;
		if(i!=frame_start){
			pos_diff = (r->p_ - pos_last)*resize0;
			pos_last = r->p_;
			rot_y_diff = rot_y-rot_y_last;
			rot_y_last = rot_y;
		}else{
			pos_diff = THR(0, 0, 0);
			rot_y_diff = 0;
			pos_last = r->p_;
			rot_y_last = rot_y;
		}
		r->p_ = THR(0, 0, 0);
		r->update_pos();
		vector<THR> d = r->get_serialized_angle();

		fprintf(f_input, "%f,%f\n", weight, weight);
		for(int j=1;j<59+1;++j){
			THR t = qua::q2em(d.at(j).q());
			fprintf(f_output, "%f,%f,%f,", t.x_, t.y_, t.z_);
		}
		fprintf(f_output, "%f,%f,%f,%f\n", pos_diff.x_, pos_diff.y_, pos_diff.z_, rot_y_diff*0.00001);
	}
	fclose(f_input);
	fclose(f_output);

	delete(r);
	return 0;
}

int make_nn_input_v()
{
	BVH bvh;
	bvh.load("./data/scene1.bvh");

	FILE* f_nn_input = NULL;
	FILE* f_nn_output = NULL;
	if((f_nn_input = fopen("f_nn_input.txt", "w")) == NULL){
		printf("file open error1.\n");
		return 1;
	}
	if((f_nn_output = fopen("f_nn_output.txt", "w")) == NULL){
		printf("file open error2.\n");
		fclose(f_nn_input);
		return 1;
	}

	vector<int> frame_feature_stick;
	//frame_feature_stick.push_back(0);
	//frame_feature_stick.push_back(1);
	frame_feature_stick.push_back(30);
	frame_feature_stick.push_back(60);
	frame_feature_stick.push_back(90);
	frame_feature_stick.push_back(120);
	frame_feature_stick.push_back(150);
	frame_feature_stick.push_back(180);
	//for(int i=0;i<10;++i){
	//	frame_feature_stick.push_back((i+1)*30);
	//}
	//frame_feature_stick.push_back(30);

	int frame_feature_max = frame_feature_stick.at(frame_feature_stick.size()-1);
	double len_stick = 0.3;

	int frame_skip = 0;
	int frame_next = 5;
	int frame_start = 0;
	int frame_end = bvh.motion_.size()/(frame_skip+1);
	//int frame_end = 241;

	vector<ROT2*> rots;
	vector<THR> pos_diff;
	vector<vector<THR> > v_part;
	vector<THR> pos_part_last;
	vector<double> y_rot_diff;
	vector<double> y_root;
	vector<pair<THR, THR> > sticks;
	THR init_pos(0, 0, 0);
	THR last_pos(0, 0, 0);
	double init_y_rot = 0;
	double last_y_rot = 0;

	int count_skip = 0;
	int i=0;
	for(auto ite = bvh.motion_.begin();ite != bvh.motion_.end();++ite, ++i, ++count_skip){
		if(count_skip >= frame_skip){
			count_skip = 0;
		}else{
			continue;
		}
		ROT2* r = ROT2::make_bone(&bvh);
		double resize0 = r->normalize_height();
		r->update_pos();
		vector<THR> data = ROT2::get_frame(bvh, i);
		r->set_serialized_angle(data);
		auto except_y = r->except_y_rotation(last_y_rot);
		double y_rot = except_y.first;
		r->q_al_cl_ = except_y.second;

		if(rots.size() == 0){
			last_pos = r->p_;
			init_pos = r->p_;
			last_y_rot = y_rot;
			init_y_rot = y_rot;
		}
		Q q = qua::e2q(-y_rot, 0, 0, qua::RotSeq::yxz);
		THR pos = (r->p_-last_pos)*resize0;
		pos = THR(q*pos.q()/q);
		pos_diff.push_back(pos);
		last_pos = r->p_;
		r->p_ = THR(0, 0, 0);
		r->update_pos();

		auto pos_part = r->get_serialized_pos();
		if(v_part.empty()){
			vector<THR> diff;
			for(int i=0;i<pos_part.size();++i){
				diff.push_back(THR(0, 0, 0));
			}
			v_part.push_back(diff);
		}else{
			vector<THR> diff;
			for(int i=0;i<pos_part.size();++i){
				diff.push_back(pos_part.at(i)-pos_part_last.at(i));
			}
			v_part.push_back(diff);
		}
		pos_part_last = pos_part;

		sticks.push_back(get_stick(r, len_stick));
		vector<double> minmax = r->get_minmax_pos();
		y_root.push_back(-minmax.at(2));
		
		y_rot_diff.push_back(y_rot-last_y_rot);
		last_y_rot = y_rot;

		rots.push_back(r);
	}

	//###nn
	for(int i=frame_start;i<frame_end-frame_feature_max;++i){
		for(int j=0;j<(int)frame_feature_stick.size();++j){
			int frame_feature = frame_feature_stick.at(j);
			THR pos_sum(0, 0, 0);
			double rot_sum = 0;
			for(int k=i;k<=i+frame_feature;++k){
				pos_sum = pos_sum+pos_diff.at(k);
				rot_sum = rot_sum+y_rot_diff.at(k);
			}
			Q q = qua::e2q(rot_sum, 0, 0, qua::RotSeq::yxz);
			THR pos_stick_f = THR(q*sticks.at(i+frame_feature).first.q()/q)+pos_sum;
			THR em_stick_f = qua::q2em(q*sticks.at(i+frame_feature).second.q());
			
			vector<THR> angle_input = rots.at(i)->get_serialized_angle();
			//vector<THR> angle_input = rots.at(i)->get_serialized_angle_al_cw();
			vector<THR> em_input;
			for(auto ite = angle_input.begin();ite != angle_input.end();++ite){
				if(ite != angle_input.begin()){
					//skip first pos_data
					em_input.push_back(qua::q2em(ite->q()));
					double y = 0;
					double x = 0;
					double z = 0;
					qua::q2e(ite->q(), y, x, z, qua::RotSeq::yxz);
					//printf("y_extract3(%f, %f, %f)\n", x*180/M_PI, y*180/M_PI, z*180/M_PI);
				}
			}
			int ii=0;
			//for unity humanoid. Needs al_cw.
			//vector<THR> angle_output = rots.at(i+frame_next)->get_serialized_angle_al_cw();
			vector<THR> angle_output = rots.at(i+frame_next)->get_serialized_angle();
			vector<THR> em_output;
			for(auto ite = angle_output.begin();ite != angle_output.end();++ite, ++ii){
				if(ite != angle_output.begin()){
					double y = 0;
					double x = 0;
					double z = 0;
					qua::q2e(ite->q(), y, x, z, qua::RotSeq::yxz);
					if(ii==1){
						printf("y_extract2(%f, %f, %f)\n", x*180/M_PI, y*180/M_PI, z*180/M_PI);
						ite->print();
					}
						
					//skip first pos_data
					THR t = *ite;
					THR em = qua::q2em(t.q());
					THR qq = qua::em2q(em);
					em_output.push_back(em);
				}
			}
			double y_root_input = y_root.at(i);
			THR pos_output = THR(0, 0, 0);
			double y_rot_output = 0;
			for(int k=i;k<=i+frame_next;++k){
				pos_output = pos_output+pos_diff.at(k);
				y_rot_output = y_rot_output+y_rot_diff.at(k);
			}
			y_rot_output = qua::convert_single_pi(y_rot_output);

			write_nn_input_v(f_nn_input, em_input, y_root_input, pos_stick_f, em_stick_f, v_part.at(i));
			write_nn_output(f_nn_output, em_output, pos_output, y_rot_output);
			printf("%d,%d,%d\n", i, j, frame_feature);
		}
	}

	for(auto ite = rots.begin();ite != rots.end();++ite){
		delete(*ite);
	}

	fclose(f_nn_input);
	fclose(f_nn_output);
	return 0;
}

int make_nn_input_traj()
{
	BVH bvh;
	bvh.load("./data/scene1.bvh");

	FILE* f_nn_input = NULL;
	FILE* f_nn_output = NULL;
	if((f_nn_input = fopen("f_nn_input.txt", "w")) == NULL){
		printf("file open error1.\n");
		return 1;
	}
	if((f_nn_output = fopen("f_nn_output.txt", "w")) == NULL){
		printf("file open error2.\n");
		fclose(f_nn_input);
		return 1;
	}

	int frame_previous = 12;
	vector<pair<THR, THR> > sticks_cw;
	vector<double> height_stick;
	double len_stick = 0.3;

	int frame_skip = 1;
	int frame_start = 0;
	int frame_end = bvh.motion_.size()/(frame_skip+1);

	vector<ROT2*> rots;

	vector<vector<THR> > v_part;
	vector<THR> pos_part_last;
	int count_skip = 0;
	int index_rots=0;
	int index_frame=0;
	for(auto ite = bvh.motion_.begin();ite != bvh.motion_.end();++ite, ++index_frame, ++count_skip){
		if(count_skip >= frame_skip){
			count_skip = 0;
		}else{
			continue;
		}
		ROT2* r = ROT2::make_bone(&bvh);
		double resize0 = r->normalize_height();
		r->update_pos();
		vector<THR> data = ROT2::get_frame(bvh, index_frame);
		r->set_serialized_angle(data);
		r->p_ = r->p_*resize0;
		r->update_pos();
		rots.push_back(r);
		
		auto st0 = get_stick(r, len_stick);
		sticks_cw.push_back(st0);

		vector<double> minmax = r->get_minmax_pos();
		height_stick.push_back(st0.first.y_ - minmax.at(2));

		if(rots.size() < frame_previous+1){
			++index_rots;
			continue;
		}

//##make_data
		ROT2* r1 = rots.at(index_rots-1);
		auto except_y = r1->except_y_rotation(rots.at(index_rots-2)->except_y_rotation().first);
		double y_rot = except_y.first;

		Q q_y_rot = qua::e2q(-y_rot, 0, 0, qua::RotSeq::yxz);
		vector<THR> pos_stick_previous;
		vector<THR> em_stick_previous;
		vector<double> height_stick_previous;
		auto st1 = sticks_cw.at(index_rots-1);
		for(int j=1;j<=frame_previous;++j){
			auto st_x = sticks_cw.at(index_rots-j);
			THR pos = st_x.first - st1.first;
			pos = THR(q_y_rot*pos.q()/q_y_rot);
			THR em = qua::q2em(q_y_rot * st_x.second.q());
			pos_stick_previous.push_back(pos);
			em_stick_previous.push_back(em);
			height_stick_previous.push_back(height_stick.at(index_rots-j));
		}

		vector<THR> v_part;
		auto pos_part1 = rots.at(index_rots-1)->get_serialized_pos();
		auto pos_part2 = rots.at(index_rots-2)->get_serialized_pos();
		for(int j=0;j<pos_part1.size();++j){
			THR p = pos_part1.at(j) - pos_part2.at(j);
			p = THR(q_y_rot*p.q()/q_y_rot);
			v_part.push_back(p);
		}

		vector<THR> angle_input = rots.at(index_rots-1)->get_serialized_angle();
		vector<THR> em_input;
		int i2 = 0;
		for(auto ite2 = angle_input.begin();ite2 != angle_input.end();++ite2, ++i2){
			if(i2 > 0){
				//skip pos
				if(i2 == 1){
					//hips
					em_input.push_back(qua::q2em(q_y_rot * ite2->q()));
				}else{
					em_input.push_back(qua::q2em(ite2->q()));
				}
			}
		}

		vector<THR> angle_output = rots.at(index_rots)->get_serialized_angle();
		vector<THR> em_output;
		i2 = 0;
		for(auto ite2 = angle_output.begin();ite2 != angle_output.end();++ite2, ++i2){
			if(i2 > 0){
				//skip pos
				if(i2 == 1){
					//hips
					em_output.push_back(qua::q2em(q_y_rot * ite2->q()));
				}else{
					em_output.push_back(qua::q2em(ite2->q()));
				}
			}
		}

		THR pos_root_diff = THR(q_y_rot*(rots.at(index_rots)->p_ - rots.at(index_rots-1)->p_).q()/q_y_rot);
		double rot_y_root_diff = rots.at(index_rots)->except_y_rotation(y_rot).first - y_rot;

//##write input
		for(auto ite2 = pos_stick_previous.begin();ite2 != pos_stick_previous.end();++ite2){
			if(ite2 != pos_stick_previous.begin()){
				fprintf(f_nn_input, "%f,%f,%f,", ite2->x_, ite2->y_, ite2->z_);
			}
		}
		for(auto ite2 = em_stick_previous.begin();ite2 != em_stick_previous.end();++ite2){
			fprintf(f_nn_input, "%f,%f,%f,", ite2->x_, ite2->y_, ite2->z_);
		}
		for(auto ite2 = height_stick_previous.begin();ite2 != height_stick_previous.end();++ite2){
/*
			if(ite2 == height_stick_previous.begin()){
				fprintf(f_nn_input, "%f", *ite2);
			}else{
				fprintf(f_nn_input, ",%f", *ite2);
			}
*/
			fprintf(f_nn_input, "%f,", *ite2);
		}

		int size_v_part_unity = 47;
		int index_unity_bone_exist[] = {
			7,9,11,36,37,38,39,40,41,42,44,45,46,48,49,50,52,53,54,56,57,58,13,14,15,16,17,18,19,21,22,23,25,26,27,29,30,31,33,34,35,4,5,6,1,2,3
		};
		for(int j=0;j<size_v_part_unity;++j){
			THR t = v_part.at(index_unity_bone_exist[j]);
			fprintf(f_nn_input, "%f,%f,%f,", t.x_, t.y_, t.z_);
		}
		for(auto ite2 = em_input.begin();ite2 != em_input.end();++ite2){
			if(ite2 == em_input.begin()){
				fprintf(f_nn_input, "%f,%f,%f", ite2->x_, ite2->y_, ite2->z_);
			}else{
				fprintf(f_nn_input, ",%f,%f,%f", ite2->x_, ite2->y_, ite2->z_);
			}
		}

		fprintf(f_nn_input, "\n");

//##write output
		{
			int ii=0;
			for(auto ite2 = em_output.begin();ite2 != em_output.end();++ite2, ++ii){
				if(ii==0){
					//skip hip
				}else if(ii==1){
					fprintf(f_nn_output, "%f,%f,%f", ite2->x_, ite2->y_, ite2->z_);
				}else{
					fprintf(f_nn_output, ",%f,%f,%f", ite2->x_, ite2->y_, ite2->z_);
				}
			}
		}
		//fprintf(f_nn_output, "%f,%f,%f,%f\n", pos_root_diff.x_, pos_root_diff.y_, pos_root_diff.z_, rot_y_root_diff);
		fprintf(f_nn_output, "\n");

		++index_rots;

	}

	for(auto ite = rots.begin();ite != rots.end();++ite){
		delete(*ite);
	}

	fclose(f_nn_input);
	fclose(f_nn_output);
	return 0;
}

int make_nn_input_one()
{
	BVH bvh;
	bvh.load("./data/scene1.bvh");

	FILE* f_nn_input = NULL;
	FILE* f_nn_output = NULL;
	if((f_nn_input = fopen("f_nn_input.txt", "w")) == NULL){
		printf("file open error1.\n");
		return 1;
	}
	if((f_nn_output = fopen("f_nn_output.txt", "w")) == NULL){
		printf("file open error2.\n");
		fclose(f_nn_input);
		return 1;
	}

	int num_dot = 100;
	double minx = -M_PI/2;
	double maxx = M_PI/2;
	double minh = -0.13;
	double maxh = 0.606;
	sm img(num_dot, num_dot);

	int frame_previous = 12;
	vector<pair<THR, THR> > sticks_cw;
	vector<double> height_stick;
	double len_stick = 0.01;

	int frame_skip = 59;
	int frame_start = 0;
	int frame_end = bvh.motion_.size()/(frame_skip+1);

	vector<ROT2*> rots;

	vector<vector<THR> > v_part;
	vector<THR> pos_part_last;
	int count_skip = 0;
	int index_rots=0;
	int index_frame=0;

	double min_height = DBL_MAX;
	double max_height = DBL_MIN;
	double min_rot_horizontal = DBL_MAX;
	double max_rot_horizontal = DBL_MIN;

	for(auto ite = bvh.motion_.begin();ite != bvh.motion_.end();++ite, ++index_frame, ++count_skip){
		if(count_skip >= frame_skip){
			count_skip = 0;
		}else{
			continue;
		}
		ROT2* r = ROT2::make_bone(&bvh);
		double resize0 = r->normalize_height();
		r->update_pos();
		vector<THR> data = ROT2::get_frame(bvh, index_frame);
		r->set_serialized_angle(data);
		r->p_ = r->p_*resize0;
		r->update_pos();
		rots.push_back(r);
		
		auto st0 = get_stick(r, len_stick);
		sticks_cw.push_back(st0);

		vector<double> minmax = r->get_minmax_pos();
		height_stick.push_back(st0.first.y_ - minmax.at(2));

		if(rots.size() < frame_previous+1){
			++index_rots;
			continue;
		}

//##make_data
		ROT2* r1 = rots.at(index_rots-1);
		auto except_y = r1->except_y_rotation(rots.at(index_rots-2)->except_y_rotation().first);
		double y_rot = except_y.first;

		Q q_y_rot = qua::e2q(-y_rot, 0, 0, qua::RotSeq::yxz);
		vector<THR> pos_stick_previous;
		vector<THR> em_stick_previous;
		vector<double> height_stick_previous;
		auto st1 = sticks_cw.at(index_rots-1);
		for(int j=1;j<=frame_previous;++j){
			auto st_x = sticks_cw.at(index_rots-j);
			THR pos = st_x.first - st1.first;
			pos = THR(q_y_rot*pos.q()/q_y_rot);
			THR em = qua::q2em(q_y_rot * st_x.second.q());
			pos_stick_previous.push_back(pos);
			em_stick_previous.push_back(em);
			height_stick_previous.push_back(height_stick.at(index_rots-j));
		}

		vector<THR> v_part;
		auto pos_part1 = rots.at(index_rots-1)->get_serialized_pos();
		auto pos_part2 = rots.at(index_rots-2)->get_serialized_pos();
		for(int j=0;j<pos_part1.size();++j){
			THR p = pos_part1.at(j) - pos_part2.at(j);
			p = THR(q_y_rot*p.q()/q_y_rot);
			v_part.push_back(p);
		}

		vector<THR> angle_input = rots.at(index_rots-1)->get_serialized_angle();
		vector<THR> em_input;
		int i2 = 0;
		for(auto ite2 = angle_input.begin();ite2 != angle_input.end();++ite2, ++i2){
			if(i2 > 0){
				//skip pos
				if(i2 == 1){
					//hips
					em_input.push_back(qua::q2em(q_y_rot * ite2->q()));
				}else{
					em_input.push_back(qua::q2em(ite2->q()));
				}
			}
		}

		vector<THR> angle_output = rots.at(index_rots)->get_serialized_angle();
		vector<THR> em_output;
		i2 = 0;
		for(auto ite2 = angle_output.begin();ite2 != angle_output.end();++ite2, ++i2){
			if(i2 > 0){
				//skip pos
				if(i2 == 1){
					//hips
					em_output.push_back(qua::q2em(q_y_rot * ite2->q()));
				}else{
					em_output.push_back(qua::q2em(ite2->q()));
				}
			}
		}

		THR pos_root_diff = THR(q_y_rot*(rots.at(index_rots)->p_ - rots.at(index_rots-1)->p_).q()/q_y_rot);
		double rot_y_root_diff = rots.at(index_rots)->except_y_rotation(y_rot).first - y_rot;

//##write input

/*
		for(auto ite2 = pos_stick_previous.begin();ite2 != pos_stick_previous.end();++ite2){
			if(ite2 != pos_stick_previous.begin()){
				fprintf(f_nn_input, "%f,%f,%f,", ite2->x_, ite2->y_, ite2->z_);
			}
		}
		for(auto ite2 = em_stick_previous.begin();ite2 != em_stick_previous.end();++ite2){
			fprintf(f_nn_input, "%f,%f,%f,", ite2->x_, ite2->y_, ite2->z_);
		}
		for(auto ite2 = height_stick_previous.begin();ite2 != height_stick_previous.end();++ite2){
			if(ite2 == height_stick_previous.begin()){
				fprintf(f_nn_input, "%f,", *ite2);
			}else{
				fprintf(f_nn_input, "%f,", *ite2);
			}
			//fprintf(f_nn_input, "%f,", *ite2);
		}
*/
/*
		int size_v_part_unity = 48;
		int index_unity_bone_exist[] = {
			0,7,9,11,36,37,38,39,40,41,42,44,45,46,48,49,50,52,53,54,56,57,58,13,14,15,16,17,18,19,21,22,23,25,26,27,29,30,31,33,34,35,4,5,6,1,2,3
		};
		for(int j=0;j<size_v_part_unity;++j){
			THR t = v_part.at(index_unity_bone_exist[j]);
			fprintf(f_nn_input, "%f,%f,%f,", t.x_, t.y_, t.z_);
		}
*/
		{
			int ii=0;
			for(auto ite2 = em_input.begin();ite2 != em_input.end();++ite2, ++ii){
				if(ii==0){
					//skip hip
				}else if(ii==1){
					fprintf(f_nn_input, "%f,%f,%f", ite2->x_, ite2->y_, ite2->z_);
				}else{
					fprintf(f_nn_input, ",%f,%f,%f", ite2->x_, ite2->y_, ite2->z_);
				}
			}
		}

		{
			//Q q = qua::em2q(em_stick_previous.at(0));
			Q q = rots.at(index_rots)->q_al_cl_.q();
			THR t = THR(q*THR(0, 0, 1).q()/q);
			double theta = asin(t.y_);
			//t.print("t");
			double height = height_stick_previous.at(0);
			fprintf(f_nn_input, ",%f,%f\n", theta, height);

			if(min_height > height) min_height = height;
			if(max_height < height) max_height = height;
			if(min_rot_horizontal > theta) min_rot_horizontal = theta;
			if(max_rot_horizontal < theta) max_rot_horizontal = theta;
		}

		//fprintf(f_nn_input, "\n");

//##write output
		{
			int ii=0;
			for(auto ite2 = em_output.begin();ite2 != em_output.end();++ite2, ++ii){
				if(ii>0){
					//skip hip
					if(ii==1){
						fprintf(f_nn_output, "%f,%f,%f", ite2->x_, ite2->y_, ite2->z_);
					}else{
						fprintf(f_nn_output, ",%f,%f,%f", ite2->x_, ite2->y_, ite2->z_);
					}
				}
			}
		}
		fprintf(f_nn_output, "\n");

		++index_rots;

	}

	printf("height_range(%f, %f)\nrot_h_range(%f, %f)\n", min_height, max_height, min_rot_horizontal, max_rot_horizontal);




	for(auto ite = rots.begin();ite != rots.end();++ite){
		delete(*ite);
	}

	fclose(f_nn_input);
	fclose(f_nn_output);
	return 0;
}
int main()
{
	srand(time(NULL));
	//t4();
	//for(int i=0;i<100;++i)
	//t8();
	//make_nn_input_one();
	make_nn_input();
	//make_nn_input_v();
	//t13();
	//t6();
	//t_rotation_order();
	//viewer();
	//t9();
	//tmp();
	//t11();
	return 0;
}
