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
	bvh.load("./data/yzx");
	//bvh.load("./data/zai");
	//bvh.load("./data/01_01.bvh");
	ROT2* r = ROT2::make_bone(&bvh);
	r->normalize_height();
	r->multiply_len(80);
	r->update_pos();
	for(int i=0;i<100;++i){
		vector<THR> data = ROT2::get_frame(bvh, i*10);
		r->set_serialized_angle(data);
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
		CamPol::simple_draw(r, buf, pos_camera, dir_camera);
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
		CamPol::simple_draw(r, buf, pos_camera, dir_camera);
	}
	delete(r);
	return 0;
}

//rots load
int t6()
{
	BVH bvh;
	//bvh.load("./data/01_01.bvh");
	bvh.load("./data/yzx");
	//bvh.load("./data/zai");
	std::vector<ROT2*> rots;
	std::vector<THR> pos_diff;
	std::vector<double> y_rotation_diff;

	THR init_pos(0, 0, 0);
	double init_y_rotation = 0;
	THR last_pos(0, 0, 0);
	double last_y_rotation = 0;

	int i=0;
	for(auto ite = bvh.motion_.begin();ite != bvh.motion_.end();++ite, ++i){
		ROT2* r = ROT2::make_bone(&bvh);
		double resize0 = r->normalize_height();
		r->update_pos();
		vector<THR> data = ROT2::get_frame(bvh, i);
		r->set_serialized_angle(data);
		if(i==0){
			init_pos = r->p_;
			init_y_rotation = r->except_y_rotation();
			printf("init_y_rotation = %f\n", init_y_rotation);
			last_pos = init_pos;
			last_y_rotation = init_y_rotation;
		}
		//r->p_ = (r->p_-init_pos)*resize0;
		last_pos = r->p_;
		//r->p_ = THR(q*r->p_.q()/q);
		//r->p_ = THR(0, 0, 0);
		double y_rotate = r->except_y_rotation();
		Q q = qua::get_quaternion_from_axis(0, 1, 0, -y_rotate);
		THR p2 = THR(q*r->p_.q()/q);
		THR init_pos2 = THR(q*init_pos.q()/q);
		r->p_ = (p2-init_pos2)*resize0;
		pos_diff.push_back(r->p_);
		y_rotate -= init_y_rotation;
		y_rotation_diff.push_back(y_rotate-last_y_rotation);
		last_y_rotation = y_rotate;

		rots.push_back(r);
	}

	for(int i=0;i<rots.size();++i){
		pos_diff.at(i).print();
		printf("%f\n", y_rotation_diff.at(i));
	}

	for(int i=0;i<100;++i){
		ROT2* r = rots.at(i*10);
		r->update_pos();
		char buf[80];
		sprintf(buf, "./out/%03d.bmp", i);
		//THR dest_camera = r->p_;
		THR dest_camera(0, 0, 0);
		THR pos_camera = THR(0, 1, 3);
		//pos_camera = THR(q*pos_camera.q()/q);
		THR dir_camera = (dest_camera-pos_camera).normalize();
		CamPol::simple_draw(r, buf, pos_camera, dir_camera);
	}

	for(auto ite = rots.begin();ite != rots.end();++ite){
		delete(*ite);
	}
	return 0;
}
int main()
{
	srand(time(NULL));
	//t4();
	//for(int i=0;i<2000;++i)
	//viewer();
	t6();
	return 0;
}
