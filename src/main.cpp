#include <stdio.h>
#include <math.h>
#include <vector>
#include <boost/tuple/tuple.hpp>
#include "bvh.h"
#include "qua.h"
#include "thr.h"
#include "rot2.h"
#include "campol.h"

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
	bvh.load("./data/a_4Char00.bvh");
	//bvh.load("./data/zai");
	//bvh.load("./data/01_01.bvh");
	ROT2* r = ROT2::make_bone(&bvh);
	r->normalize_height();
	r->multiply_len(80);
	r->update_pos();
	double rot_init = 0;
	for(int i=0;i<100;++i){
		vector<THR> data = ROT2::get_frame(bvh, i*10);
		r->set_serialized_angle(data);
		double rot = r->except_y_rotation();
		if(i==0) rot_init = rot;
		Q qq = qua::e2q(rot, 0, 0, qua::RotSeq::yxz);
		r->q_al_cl_ = THR(qq*r->q_al_cl_.q());
		printf("rot_=%f\n", rot-rot_init);
		r->update_pos();
		//r->print();
		char buf[80];
		sprintf(buf, "./out/%03d.bmp", i);
		//THR dest_camera = r->p_;
		THR dest_camera(0, 100, 0);
		Q q = qua::e2q(0, i*M_PI/50, 0, qua::RotSeq::xyz);
		THR pos_camera = THR(160, 120, 0);
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
		//r->p_ = THR(0, 0, 0);
		double y_rot = 0;
		//r->p_.print();
		y_rot = r->except_y_rotation();
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
	bvh.load("./data/a_4Char00.bvh");

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
		double y_rot = 0;
		y_rot = r->except_y_rotation();
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

	for(int i=0;i<rots.size();++i){
		//pos_diff.at(i).print();
		//printf("%f\n", y_rotation_diff.at(i));
		//printf("y_root = %f\n", y_root.at(i));
	}

	THR pos_sum(0, 0, 0);
	double rot_sum = 0;//init_y_rot;
	for(int i=0;i<100;++i){
		for(int j=0;j<10;++j){
			if(i == 0) break;
			pos_sum = pos_sum + pos_diff.at((i-1)*10+j);
			rot_sum = rot_sum + y_rot_diff.at((i-1)*10+j);
		}
		ROT2* r = rots.at(i*10);
		r->p_ = r->p_+pos_sum;
		Q q = qua::e2q(rot_sum, 0, 0, qua::RotSeq::yxz);
		r->q_al_cl_ = THR(q*r->q_al_cl_.q());
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
	double rot_init = 0;
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
	bvh.load("./data/a_4Char00.bvh");

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
	frame_feature_stick.push_back(60);
	frame_feature_stick.push_back(120);
	frame_feature_stick.push_back(240);
	int frame_feature_max = frame_feature_stick.at(frame_feature_stick.size()-1);
	double len_stick = 0.3;

	vector<ROT2*> rots;
	vector<THR> pos_diff;
	vector<double> y_rot_diff;
	vector<double> y_root;
	vector<pair<THR, THR> > sticks;
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
		double y_rot = 0;
		y_rot = r->except_y_rotation();
		if(i==0){
			last_pos = r->p_;
			init_pos = r->p_;
			last_y_rot = y_rot;
			init_y_rot = y_rot;
		}
		Q q = qua::e2q(-init_y_rot, 0, 0, qua::RotSeq::yxz);
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

	//nn
	int frame_next = 10;
	int frame_start = 0;
	int frame_end = rots.size();
	for(int i=frame_start;i<frame_end-frame_feature_max;++i){
		for(int j=0;j<frame_feature_stick.size();++j){
			int frame_feature = frame_feature_stick.at(j);
			THR pos_sum(0, 0, 0);
			double rot_sum = 0;
			for(int k=i;k<i+frame_feature;++k){
				pos_sum = pos_sum+pos_diff.at(k);
				rot_sum = rot_sum+y_rot_diff.at(k);
			}
			Q q = qua::e2q(rot_sum, 0, 0, qua::RotSeq::yxz);
			THR pos_stick_f = THR(q*sticks.at(i+frame_feature).first.q()/q)+pos_sum;
			THR em_stick_f = qua::q2em(q*sticks.at(i+frame_feature).second.q());
			
			vector<THR> angle_input = rots.at(i)->get_serialized_angle();
			vector<THR> em_input;
			for(auto ite = angle_input.begin();ite != angle_input.end();++ite){
				if(ite != angle_input.begin()){
					//skip first pos_data
					em_input.push_back(qua::q2em(ite->q()));
				}
			}
			vector<THR> angle_output = rots.at(i+frame_next)->get_serialized_angle();
			vector<THR> em_output;
			for(auto ite = angle_output.begin();ite != angle_output.end();++ite){
				if(ite != angle_output.begin()){
					//skip first pos_data
					em_output.push_back(qua::q2em(ite->q()));
				}
			}
			double y_root_input = y_root.at(i);
			THR pos_output = THR(0, 0, 0);
			double y_rot_output = 0;
			for(int k=i;k<i+frame_next;++k){
				pos_output = pos_output+pos_diff.at(k);
				y_rot_output = y_rot_output+y_rot_diff.at(k);
			}

			write_nn_input(f_nn_input, em_input, y_root_input, pos_stick_f, em_stick_f);
			write_nn_output(f_nn_output, em_output, pos_output, y_rot_output);
			printf("%d,%d,%d\n", i, j, frame_feature);
		}
	}

	//draw
	if(false){
	THR pos_sum(0, 0, 0);
	double rot_sum = 0;//init_y_rot;
	for(int i=0;i<100;++i){
		for(int j=0;j<10;++j){
			if(i == 0) break;
			pos_sum = pos_sum + pos_diff.at((i-1)*10+j);
			rot_sum = rot_sum + y_rot_diff.at((i-1)*10+j);
		}
		ROT2* r = rots.at(i*10);
		r->p_ = r->p_+pos_sum;
		Q q = qua::e2q(rot_sum, 0, 0, qua::RotSeq::yxz);
		r->q_al_cl_ = THR(q*r->q_al_cl_.q());
		r->update_pos();
		char buf[80];
		sprintf(buf, "./out/%03d.bmp", i);
		THR dest_camera(0, 0, 0);
		THR pos_camera = THR(0, 1, 3);
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
/*
int make_nn_input()
{
	vector<int> frame_feature_stick;
	frame_feature_stick.push_back(60);
	frame_feature_stick.push_back(120);
	frame_feature_stick.push_back(240);
	int frame_feature_max = frame_feature_stick.at(frame_feature_stick.size()-1);
	double len_stick = 0.3;

	BVH bvh;
	bvh.load("./data/scene1.bvh");

	ROT2* r = ROT2::make_bone(&bvh);
	double resize0 = r->normalize_height();
	r->update_pos();
	double rot_init = 0;
	for(int i=0;i<100;++i){
		vector<THR> data = ROT2::get_frame(bvh, i*10);
		r->set_serialized_angle(data);
		r->p_ = r->p_*resize0*resize1;
		r->update_pos();
		//r->print();
		char buf[80];
		sprintf(buf, "./out/%03d.bmp", i);
		//THR dest_camera = r->p_;
		THR dest_camera(0, 100, 0);
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
*/
int main()
{
	srand(time(NULL));
	//t4();
	//for(int i=0;i<100;++i)
	//t8();
	make_nn_input();
	//viewer();
	return 0;
}
