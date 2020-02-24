/*
 
     FreeSurround Output Plugin
      
    Copyright (c) 2009 Michel Cailhol
    Based on foo_dsp_freesurround (c) 2007 Christian Kothe
     
    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */

#include <stdio.h>
#include <string.h>
#define __USE_XOPEN
#include <unistd.h>
#include <alsa/asoundlib.h>
#include <alsa/pcm_external.h>
#include <alsa/pcm_plugin.h>
#include <fftw3.h>
#include <complex.h>
#include <math.h>

const float PI = 3.141592654;
const float epsilon = 0.000001;
const float center_level = 0.353553385; // 0.5*sqrt(0.5);	// gain of center channel
const int BLOCK_SIZE = 8192;
const unsigned int DBL_BLOCK_SIZE = 2 * 8192;
const unsigned int DBL_BLOCK_SIZE_MIN_ONE = (2 * 8192) - 1;
const int INPUT_CHANNELS = 2;
const int OUTPUT_CHANNELS = 6;



typedef float complex cfloat;
typedef float farraybs[8192];
typedef float farraybse[12288]; // BLOCK_SIZE + BLOCK_SIZE/2
typedef short sarraydbs[16384];
typedef cfloat cfarraybs[8192];

typedef struct snd_pcm_fsupmix snd_pcm_fsupmix_t;

typedef void (*fsupmixer_t)(snd_pcm_fsupmix_t *fsmix,
			  const snd_pcm_channel_area_t *dst_areas,
			  snd_pcm_uframes_t dst_offset,
			  const snd_pcm_channel_area_t *src_areas,
			  snd_pcm_uframes_t src_offset,
			  snd_pcm_uframes_t size);


struct snd_pcm_fsupmix {
    snd_pcm_extplug_t ext;
	snd_pcm_format_t format;
	unsigned int channels;
	unsigned int rate;
	unsigned int bitrate;
	int outbuf_size;
	snd_pcm_uframes_t transfer;
	int remain;
	unsigned int slave_period_size;
	unsigned int slave_buffer_size;
	snd_pcm_hw_params_t *hw_params;
    farraybs myinput[2];
    sarraydbs myoutput[6];
    unsigned int myout_count;
    unsigned int myout_write_offset;
    unsigned int myout_read_offset;
	char debug;

    int preupmixed;
    	// FFTW data structures
	float *lt,*rt,*dst;				   // left total, right total (source arrays), destination array
	fftwf_complex *dftL,*dftR,*src;    // intermediate arrays (FFTs of lt & rt, processing source)
	fftwf_plan loadL,loadR,store;      // plans for loading the data into the intermediate format and back

	// buffers
    cfarraybs frontL; // the signal (phase-corrected) in the frequency domain
    cfarraybs frontR; // the signal (phase-corrected) in the frequency domain
    cfarraybs avg; // the signal (phase-corrected) in the frequency domain
    cfarraybs surL; // the signal (phase-corrected) in the frequency domain
    cfarraybs surR; // the signal (phase-corrected) in the frequency domain
    cfarraybs trueavg;       // for lfe generation
    
    farraybs xfs ; // the feature space positions for each frequency bin
    farraybs yfs ; // the feature space positions for each frequency bin
    farraybs wnd ; // the window function, precalculated
	farraybs filter[6] ; // a frequency filter for each output channel    
    
    farraybse fsinbuf[2];	   // the sliding input buffers
	farraybse fsoutbuf[6];	   // the sliding output buffers  
    
        	// coefficients
	float surround_high,surround_low;  // high and low surround mixing coefficient (e.g. 0.8165/0.5774)
	float surround_balance;			   // the xfs balance that follows from the coeffs
	float surround_level;			   // gain for the surround channels (follows from the coeffs
	float master_gain;				   // gain for all channels
	float phase_offsetL, phase_offsetR;// phase shifts to be applied to the rear channels
	float front_separation;			   // front stereo separation
	float rear_separation;			   // rear stereo separation
	int linear_steering;			   // whether the steering should be linear or not 
    float center_width;                // distribution of the center information towards the front left/right channels
    float dimension;
    float adaption_rate;
    unsigned int phase_mode;
    
};



static inline void *area_addr(const snd_pcm_channel_area_t *area, snd_pcm_uframes_t offset) {
	unsigned int bitofs = area->first + area->step * offset;
	return (char *) area->addr + bitofs / 8;
}

static inline unsigned int area_step(const snd_pcm_channel_area_t *area) {
	return area->step / 8;
}




   	// set the assumed surround mixing coefficients
static	void surround_coefficients(snd_pcm_fsupmix_t *fsmix, float a, float b) {
		fsmix->master_gain = 1.0;
		// calc the simple coefficients
		fsmix->surround_high = a;
		fsmix->surround_low = b;
		fsmix->surround_balance = (a-b)/(a+b);
		// increase the rear volume
		//fsmix->surround_level = 1/(a+b);
		fsmix->surround_level = 0.85;
}


// set the phase shifting mode for decoding
	// 0 = (+0°,+0°)   - music mode
	// 1 = (+0°,+180°) - PowerDVD compatibility
	// 2 = (+180°,+0°) - BeSweet compatibility
	// 3 = (-90°,+90°) - experimental exact reconstuction; not clear whether this is entirely correct
static	void set_phase_mode(snd_pcm_fsupmix_t *fsmix) {
		const float modes[4][2] = {{0,0},{0,PI},{PI,0},{-PI/2,PI/2}};
		fsmix->phase_offsetL = modes[fsmix->phase_mode][0];
		fsmix->phase_offsetR = modes[fsmix->phase_mode][1];
}

	// what steering mode should be chosen
/*static	void steering_mode(snd_pcm_fsupmix_t *fsmix, int mode) { 
    fsmix->linear_steering = mode; 
}

	// set front & rear separation controls
static	void separation(snd_pcm_fsupmix_t *fsmix, float front, float rear) {
		fsmix->front_separation = front;
		fsmix->rear_separation = rear;
}
*/
   // set lfe filter params
static void sample_rate(snd_pcm_fsupmix_t *fsmix, unsigned int srate) {
        // lfe filter is just straight through band limited
        unsigned int cutoff = (250*BLOCK_SIZE)/srate;
        unsigned f;
        for (f=0;f<=BLOCK_SIZE/2;f++) {           
            if ((f>=2) && (f<cutoff))
                //fsmix->filter[5][f] = 0.5*sqrt(0.5);
                fsmix->filter[5][f] = 1.0;
            else
                fsmix->filter[5][f] = 0.0;
        }
}




// polar <-> cartesian coodinates conversion
static inline float amplitude(const float cf[2]) { return sqrt(cf[0]*cf[0] + cf[1]*cf[1]); }
static inline float phase(const float cf[2]) { return atan2(cf[1],cf[0]); }
static  cfloat polar(float a, float p) {
    cfloat cpf = a*cos(p) + I*(a*sin(p));
    return cpf;
     //return cfloat (a*cos(p),a*sin(p)); 
}
static inline float sqr(float x) { return x*x; }
// the dreaded min/max
static inline float min(float a, float b) { return a<b?a:b; }
static inline float max(float a, float b) { return a>b?a:b; }
static inline float clamp(float x) { return max(-1,min(1,x)); }



	// map from amplitude difference and phase difference to yfs
static	inline double get_yfs(double ampDiff, double phaseDiff) {
		double x = 1-(((1-sqr(ampDiff))*phaseDiff)/PI*2);
		return 0.16468622925824683 + 0.5009268347818189*x - 0.06462757726992101*x*x
			+ 0.09170680403453149*x*x*x + 0.2617754892323973*tan(x) - 0.04180413533856156*sqr(tan(x));
}

	// map from amplitude difference and yfs to xfs
static	inline double get_xfs(double ampDiff, double yfs) {
		double x=ampDiff,y=yfs;
		return 2.464833559224702*x - 423.52131153259404*x*y + 
			67.8557858606918*x*x*x*y + 788.2429425544392*x*y*y - 
			79.97650354902909*x*x*x*y*y - 513.8966153850349*x*y*y*y + 
			35.68117670186306*x*x*x*y*y*y + 13867.406173420834*y*asin(x) - 
			2075.8237075786396*y*y*asin(x) - 908.2722068360281*y*y*y*asin(x) - 
			12934.654772878019*asin(x)*sin(y) - 13216.736529661162*y*tan(x) + 
			1288.6463247741938*y*y*tan(x) + 1384.372969378453*y*y*y*tan(x) + 
			12699.231471126128*sin(y)*tan(x) + 95.37131275594336*sin(x)*tan(y) - 
			91.21223198407546*tan(x)*tan(y);
}

	// filter the complex source signal and add it to target
static	void apply_filter(snd_pcm_fsupmix_t *fsmix, cfarraybs signal, int nfilter, int nchannel ) {
		// filter the signal
        unsigned f;
		for (f=0; f<=BLOCK_SIZE/2; f++) {		
			fsmix->src[f][0] = crealf(signal[f]) * fsmix->filter[nfilter][f];
			fsmix->src[f][1] = cimagf(signal[f]) * fsmix->filter[nfilter][f];
		}
		// transform into time domain
		fftwf_execute(fsmix->store);
		// add the result to target, windowed
        unsigned k;
		for ( k=0; k<BLOCK_SIZE; k++) {
			fsmix->fsoutbuf[nchannel][k + BLOCK_SIZE/2] += fsmix->wnd[k]*fsmix->dst[k];
        }
}




	// CORE FUNCTION: decode a block of data
static	void block_decode(snd_pcm_fsupmix_t *fsmix, int buffoffset) {
    
    
    // 1. scale the input by the window function; this serves a dual purpose:
    // - first it improves the FFT resolution b/c boundary discontinuities (and their frequencies) get removed
    // - second it allows for smooth blending of varying filters between the blocks
    unsigned k;
    for (k=0; k<BLOCK_SIZE; k++) {
        fsmix->lt[k] = fsmix->fsinbuf[0][k + buffoffset] * fsmix->wnd[k] * fsmix->master_gain;
        fsmix->rt[k] = fsmix->fsinbuf[1][k + buffoffset] * fsmix->wnd[k] * fsmix->master_gain;
    }

    // ... and tranform it into the frequency domain
    fftwf_execute(fsmix->loadL);
    fftwf_execute(fsmix->loadR);

    // 2. compare amplitude and phase of each DFT bin and produce the X/Y coordinates in the sound field
    unsigned f;
    for (f=0; f<=BLOCK_SIZE/2; f++) {			
        // get left/right amplitudes/phases
        float ampL = amplitude(fsmix->dftL[f]), ampR = amplitude(fsmix->dftR[f]);
        float phaseL = phase(fsmix->dftL[f]), phaseR = phase(fsmix->dftR[f]);

        // calculate the amplitude/phase difference
        float ampDiff = clamp((ampL+ampR < epsilon) ? 0 : (ampR-ampL) / (ampR+ampL));
        float phaseDiff = phaseL - phaseR;
        if (phaseDiff < -PI) phaseDiff += 2*PI;
        if (phaseDiff > PI) phaseDiff -= 2*PI;
        phaseDiff = fabs(phaseDiff);

        if (fsmix->linear_steering) {
            // --- the new linear mode ---

            // get sound field x/y position
            fsmix->yfs[f] = get_yfs(ampDiff,phaseDiff);
            fsmix->xfs[f] = get_xfs(ampDiff, fsmix->yfs[f]);

            // add dimension control
            fsmix->yfs[f] = clamp(fsmix->yfs[f] - fsmix->dimension);

            // add crossfeed control
            fsmix->xfs[f] = clamp(fsmix->xfs[f] * (fsmix->front_separation*(1+fsmix->yfs[f])/2 + fsmix->rear_separation*(1-fsmix->yfs[f])/2));

            // 3. generate frequency filters for each output channel
            float left = (1-fsmix->xfs[f])/2, right = (1+fsmix->xfs[f])/2;
            float front = (1+fsmix->yfs[f])/2, back = (1-fsmix->yfs[f])/2;
            float volume[5] = {
                front * (left * fsmix->center_width + max(0,-fsmix->xfs[f]) * (1-fsmix->center_width)),	// left
                front * center_level*((1-fabs(fsmix->xfs[f])) * (1-fsmix->center_width)),			// center
                front * (right * fsmix->center_width + max(0, fsmix->xfs[f]) * (1-fsmix->center_width)),	// right
                back * fsmix->surround_level * left,										// left surround
                back * fsmix->surround_level * right										// right surround
            };

            // adapt the prior filter
            unsigned c;
            for (c=0;c<5;c++)
                fsmix->filter[c][f] = (1-fsmix->adaption_rate)*fsmix->filter[c][f] + fsmix->adaption_rate*volume[c]/BLOCK_SIZE;

        } else {
            // --- the old & simple steering mode ---
            // calculate the amplitude/phase difference
				float ampDiff = clamp((ampL+ampR < epsilon) ? 0 : (ampR-ampL) / (ampR+ampL));
				float phaseDiff = phaseL - phaseR;
				if (phaseDiff < -PI) phaseDiff += 2*PI;
				if (phaseDiff > PI) phaseDiff -= 2*PI;
				phaseDiff = fabs(phaseDiff);

				// determine sound field x-position
				fsmix->xfs[f] = ampDiff;

				// determine preliminary sound field y-position from phase difference
				fsmix->yfs[f] = 1 - (phaseDiff/PI)*2;

				if (fabs(fsmix->xfs[f]) > fsmix->surround_balance) {
					// blend linearly between the surrounds and the fronts if the balance exceeds the surround encoding balance
					// this is necessary because the sound field is trapezoidal and will be stretched behind the listener
					float frontness = (fabs(fsmix->xfs[f]) - fsmix->surround_balance)/(1-fsmix->surround_balance);
					fsmix->yfs[f]  = (1-frontness) * fsmix->yfs[f] + frontness * 1; 
				}

				// add dimension control
				fsmix->yfs[f] = clamp(fsmix->yfs[f] - fsmix->dimension);

				// add crossfeed control
				fsmix->xfs[f] = clamp(fsmix->xfs[f] * (fsmix->front_separation*(1+fsmix->yfs[f])/2 + fsmix->rear_separation*(1-fsmix->yfs[f])/2));

				// 3. generate frequency filters for each output channel, according to the signal position
				// the sum of all channel volumes must be 1.0
				float left = (1-fsmix->xfs[f])/2, right = (1+fsmix->xfs[f])/2;
				float front = (1+fsmix->yfs[f])/2, back = (1-fsmix->yfs[f])/2;
				float volume[5] = {
					front * (left * fsmix->center_width + max(0,-fsmix->xfs[f]) * (1-fsmix->center_width)),		// left
					front * center_level*((1-fabs(fsmix->xfs[f])) * (1-fsmix->center_width)),				// center
					front * (right * fsmix->center_width + max(0, fsmix->xfs[f]) * (1-fsmix->center_width)),		// right
					back * fsmix->surround_level*max(0,min(1,((1-(fsmix->xfs[f]/fsmix->surround_balance))/2))),	// left surround
					back * fsmix->surround_level*max(0,min(1,((1+(fsmix->xfs[f]/fsmix->surround_balance))/2)))	// right surround
				};

				// adapt the prior filter
                unsigned c;
				for (c=0;c<5;c++)
					fsmix->filter[c][f] = (1-fsmix->adaption_rate)*fsmix->filter[c][f] + fsmix->adaption_rate*volume[c]/BLOCK_SIZE;
            }

        // ... and build the signal which we want to position
        fsmix->frontL[f] = polar(ampL+ampR,phaseL);
        fsmix->frontR[f] = polar(ampL+ampR,phaseR);
        fsmix->avg[f] = fsmix->frontL[f] + fsmix->frontR[f];
        fsmix->surL[f] = polar(ampL+ampR,phaseL+fsmix->phase_offsetL);
        fsmix->surR[f] = polar(ampL+ampR,phaseR+fsmix->phase_offsetR);
        fsmix->trueavg[f] = (fsmix->dftL[f][0] + fsmix->dftR[f][0]) + I*(fsmix->dftL[f][1] + fsmix->dftR[f][1]);
    }

    // 4. distribute the unfiltered reference signals over the channels
    apply_filter(fsmix, fsmix->frontL, 0, 0);	// front left
    apply_filter(fsmix, fsmix->avg, 1, 4);		// front center
    apply_filter(fsmix, fsmix->frontR, 2, 1);	// front right
    apply_filter(fsmix, fsmix->surL, 3, 2);		// surround left
    apply_filter(fsmix, fsmix->surR, 4, 3);		// surround right

    /*
	just to test on stereo card
	apply_filter(fsmix, fsmix->avg, 3, 0);		// left
    apply_filter(fsmix, fsmix->avg, 4, 1);		// right*/

    apply_filter(fsmix, fsmix->trueavg, 5, 5);  // lfe
}




	// handle the output buffering for overlapped calls of block_decode
static	void add_output(snd_pcm_fsupmix_t *fsmix, int buffoffset, int result) {
    
    // add the windowed data to the last 2/3 of the output buffer
    block_decode(fsmix, buffoffset);
    unsigned c, k;
    unsigned int ooffset;
  
    for (c=0;c<6;c++) {
        if (result) 
            // return the first 2/3 of the ouput buffer
            for (k=0;k<BLOCK_SIZE;k++) {
                ooffset = (fsmix->myout_write_offset + k) & DBL_BLOCK_SIZE_MIN_ONE;
                if (c != 5)
                    fsmix->myoutput[c][ooffset] = (short)(fsmix->fsoutbuf[c][k] + 0.5);
                else // lfe channel
                    fsmix->myoutput[5][ooffset] = (fsmix->myoutput[0][ooffset] >> 1) + (fsmix->myoutput[1][ooffset] >> 1);
            }

        for ( k=0;k<BLOCK_SIZE;k++)
            // shift the last 2/3 to the first 2/3 of the output buffer
            fsmix->fsoutbuf[c][k] = fsmix->fsoutbuf[c][k+BLOCK_SIZE/2];
        // and clear the rest
        for (k=BLOCK_SIZE; k<BLOCK_SIZE+BLOCK_SIZE/2; k++)
            fsmix->fsoutbuf[c][k] = 0;
    }
    if (result) {
        fsmix->myout_count += BLOCK_SIZE;
        fsmix->myout_write_offset += BLOCK_SIZE;
        fsmix->myout_write_offset &= DBL_BLOCK_SIZE_MIN_ONE;
    }        
}




	// decode a chunk of stereo sound, has to contain exactly blocksize samples
static	void decode(snd_pcm_fsupmix_t *fsmix) {
    // append incoming data to the end of the input buffer 
    unsigned k;
    for (k=0; k<BLOCK_SIZE; k++) {		
        fsmix->fsinbuf[0][k+BLOCK_SIZE/2] = fsmix->myinput[0][k];
        fsmix->fsinbuf[1][k+BLOCK_SIZE/2] = fsmix->myinput[1][k];
    }
    // process first part
    add_output(fsmix, 0, 0);
    // process second part (overlapped) and return result
    add_output(fsmix, BLOCK_SIZE/2, 1);
    // shift last third of input buffer to the beginning
    for (k=0; k<BLOCK_SIZE/2; k++) {		
        fsmix->fsinbuf[0][k] = fsmix->fsinbuf[0][k+BLOCK_SIZE];
        fsmix->fsinbuf[1][k] = fsmix->fsinbuf[1][k+BLOCK_SIZE];
    }
}


/*
 * transfer callback
 */
static snd_pcm_sframes_t fs_transfer(snd_pcm_extplug_t *ext,
	       const snd_pcm_channel_area_t *dst_areas,
	       snd_pcm_uframes_t dst_offset,
	       const snd_pcm_channel_area_t *src_areas,
	       snd_pcm_uframes_t src_offset,
	       snd_pcm_uframes_t size)
{
	snd_pcm_fsupmix_t *fsmix = (snd_pcm_fsupmix_t *)ext;




	if (ext->format==ext->slave_format &&
		ext->subformat==ext->slave_subformat &&
		ext->channels==1)
	{

		snd_pcm_areas_copy(dst_areas+4, dst_offset, src_areas, src_offset,
				   1, size, ext->format);
		return size;
	}


    unsigned int len2 = BLOCK_SIZE - fsmix->preupmixed;
	short *src[6], *dst[6];
	unsigned int src_step[6], dst_step[6], i, k;
	

	char real_surround=0;



	if (size > len2) {
		size = len2;
       // printf("on en a trop : size2=%u, len2=%u\n",size2,len2);
    }
    //printf("size1=%u\n",size);    
    for (i = 0; i < ext->channels; i++) {
		src[i] = (short *)area_addr(src_areas + i, src_offset);
		src_step[i] = area_step(src_areas + i) / 2;
	}
    for (i = 0; i < OUTPUT_CHANNELS; i++) {
		dst[i] = (short *)area_addr(dst_areas + i, dst_offset);
    //    printf("dst[%u]=%p\n",i,dst[i]);
		dst_step[i] = area_step(dst_areas + i) / 2;
	}

	if (ext->channels>=4)
	{
		for (k=0; k<size && !real_surround; k++) 
		{
			real_surround=real_surround || *src[2] || *src[3]|| *src[4]|| *src[5];
			src[2] += src_step[2];
			src[3] += src_step[3];
			src[4] += src_step[4];
			src[5] += src_step[5];
		}
		if (real_surround)
		{
			snd_pcm_areas_copy(dst_areas, dst_offset, src_areas, src_offset,
					ext->channels, size, ext->format);
			return size;
		}
	}




	for (k=0; k<size; k++) {
        fsmix->myinput[0][fsmix->preupmixed + k] = (float)*src[0];
        fsmix->myinput[1][fsmix->preupmixed + k] = (float)*src[1];
		src[0] += src_step[0];
		src[1] += src_step[1];
	}



/*	if (ext->format==ext->slave_format &&
		ext->subformat==ext->slave_subformat &&
		ext->channels==6)
	{
	}

	if (ext->format==ext->slave_format &&
		ext->subformat==ext->slave_subformat &&
		ext->channels==ext->slave_channels)
	{
		snd_pcm_areas_copy(dst_areas, dst_offset, src_areas, src_offset,
				   ext->channels, size, ext->format);
		return size;
	}

*/
    fsmix->preupmixed += size;
    
    if (fsmix->preupmixed == BLOCK_SIZE && fsmix->myout_count > BLOCK_SIZE) {
        fsmix->preupmixed -= size;
        size = 0;
    }
    else
    if (fsmix->preupmixed == BLOCK_SIZE) { 
    //    printf("on a rempli un bloc : on ecrit a partir de %u\n",fsmix->myout_write_offset);
        decode(fsmix);  //  adaption_rate 
    /*    unsigned int ooffset;
        for (k=0; k<BLOCK_SIZE; k++) {
            ooffset = (rec->myout_write_offset + k) & DBL_BLOCK_SIZE_MIN_ONE;
            rec->myoutput[0][ooffset] = (short)(rec->myinput[0][k] + 0.5);
            rec->myoutput[1][ooffset] = (short)(rec->myinput[1][k] + 0.5);
            rec->myoutput[2][ooffset] = (short)(rec->myinput[0][k] + 0.5);
            rec->myoutput[3][ooffset] = (short)(rec->myinput[1][k] + 0.5);
            rec->myoutput[4][ooffset] = 0;
            rec->myoutput[5][ooffset] = 0;
        }
        rec->myout_count += BLOCK_SIZE;
        rec->myout_write_offset += BLOCK_SIZE;
        rec->myout_write_offset &= DBL_BLOCK_SIZE_MIN_ONE; */
        fsmix->preupmixed = 0;
    }
  //  printf("buffer contient %i\n", fsmix->myout_count);
    
    unsigned int ooffset;
    /* flatten copy to n-channel interleaved */
	double levels[6]={0,0,0,0,0,0};
    for (k = 0; k < size; k++) {
        for (i = 0; i < OUTPUT_CHANNELS; i++) {
        //    printf("k=%u, i=%u\n",k,i);
            ooffset = (fsmix->myout_read_offset + k) & DBL_BLOCK_SIZE_MIN_ONE;
            *dst[i] = fsmix->myoutput[i][ooffset];
            dst[i] += dst_step[i];   
			levels[i]=levels[i]+abs(*dst[i]);
        }
    }
    
    fsmix->myout_count -= size;
    fsmix->myout_read_offset += size;
    fsmix->myout_read_offset &= DBL_BLOCK_SIZE_MIN_ONE;


	if (fsmix->debug)
	{
		for(i=0;i<6;i++)
			for(k=0;k<20;k++)
				printf(levels[i]/size/500>k?"#":" ");
		printf("\n");
	}
	return size;
}




/*
 * prepare callback
 *
 * Allocate internal buffers
 */
static int fs_prepare(snd_pcm_extplug_t *ext)
{
	snd_pcm_fsupmix_t *fsmix = (snd_pcm_fsupmix_t *)ext;

        
    fsmix->myout_count = BLOCK_SIZE;
    fsmix->myout_write_offset = BLOCK_SIZE;
    fsmix->myout_read_offset = 0;
    
       // create FFTW buffers
    fsmix->lt = (float*)fftwf_malloc(sizeof(float)*BLOCK_SIZE);
    fsmix->rt = (float*)fftwf_malloc(sizeof(float)*BLOCK_SIZE);
    fsmix->dst = (float*)fftwf_malloc(sizeof(float)*BLOCK_SIZE);
    fsmix->dftL = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*BLOCK_SIZE);
    fsmix->dftR = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*BLOCK_SIZE);
    fsmix->src = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*BLOCK_SIZE);
    fsmix->loadL = fftwf_plan_dft_r2c_1d(BLOCK_SIZE, fsmix->lt, fsmix->dftL,FFTW_MEASURE);
    fsmix->loadR = fftwf_plan_dft_r2c_1d(BLOCK_SIZE, fsmix->rt, fsmix->dftR,FFTW_MEASURE);
    fsmix->store = fftwf_plan_dft_c2r_1d(BLOCK_SIZE, fsmix->src, fsmix->dst,FFTW_MEASURE);    
    

	fsmix->transfer = 0;
	fsmix->remain = 0;
    fsmix->preupmixed = 0;
    
    unsigned k;
    for (k=0;k<BLOCK_SIZE;k++)
        fsmix->wnd[k] = sqrt(0.5*(1-cos(2*PI*k/BLOCK_SIZE)));  // square root of hann
    
    surround_coefficients(fsmix, 0.8165,0.5774);
    set_phase_mode(fsmix);    
    //separation(fsmix, 1.0,1.0);
	//steering_mode(rec, 1);
    sample_rate(fsmix, ext->rate);

 //   printf("prepare :  center_width=%f dimension=%f adaption_rate=%f phase_mode=%u linear_steering=%i\n", 
 //           rec->center_width, rec->dimension, rec->adaption_rate, rec->phase_mode, rec->linear_steering);

	return 0;
}



/*
 * close callback
 */
static int fs_close(snd_pcm_extplug_t *ext)
{
	snd_pcm_fsupmix_t *fsmix = (snd_pcm_fsupmix_t *)ext;

        // clean up the FFTW stuff
    fftwf_destroy_plan(fsmix->store);
    fftwf_destroy_plan(fsmix->loadR);
    fftwf_destroy_plan(fsmix->loadL);
    fftwf_free(fsmix->src); 
    fftwf_free(fsmix->dftR);
    fftwf_free(fsmix->dftL);
    fftwf_free(fsmix->dst);
    fftwf_free(fsmix->rt);
    fftwf_free(fsmix->lt);;
	return 0;
}
			      
/*
 * callback table
 */
static snd_pcm_extplug_callback_t fs_callback = {
	.transfer = fs_transfer,
	.close = fs_close,
	.init = fs_prepare,
};



/*
 * Main entry point
 */
SND_PCM_PLUGIN_DEFINE_FUNC(freesurround)
{
	snd_config_iterator_t i, next;
    snd_pcm_fsupmix_t *fsmix;
    snd_config_t *sconf = NULL;
    static const unsigned int chlist[2] = {4, 6};    
	int err;
	unsigned int rate = 48000;
	unsigned int bitrate = 448;
	unsigned int channels = 6;
	snd_pcm_format_t format = SND_PCM_FORMAT_S16_LE;
    float center_width = 0.50;
    float dimension = 0.25;
    float adaption_rate = 0.9;
    unsigned int phase_mode = 0;
    int linear_steering = 1;
    float front_separation = 1.0;
    float rear_separation = 1.0;
	char debug=0;


	if (stream != SND_PCM_STREAM_PLAYBACK) {
		SNDERR("freesurround is only for playback");
		return -EINVAL;
	}

	snd_config_for_each(i, next, conf) {
		snd_config_t *n = snd_config_iterator_entry(i);
		const char *id;
		if (snd_config_get_id(n, &id) < 0)
			continue;
		if (strcmp(id, "comment") == 0 || strcmp(id, "type") == 0 || strcmp(id, "hint") == 0)
			continue;
		if (strcmp(id, "slave") == 0) {
			sconf = n;
			continue;
		}

		if (strcmp(id, "channels") == 0) {
			long val;
			if (snd_config_get_integer(n, &val) < 0) {
				SNDERR("Invalid type for %s", id);
				return -EINVAL;
			}
			channels = val;
			if (channels != 2 && channels != 4 && channels != 6) {
				SNDERR("channels must be 2, 4 or 6");
				return -EINVAL;
			}
			continue;
		}

		if (strcmp(id, "debug") == 0) {
			long val;
			if (snd_config_get_integer(n, &val) < 0) {
				SNDERR("Invalid type for %s", id);
				return -EINVAL;
			}
			debug = val;
			if (debug != 0 && debug != 1 ) {
				SNDERR("debug must be 0 or 1");
				return -EINVAL;
			}
			continue;
		}


        if (strcmp(id, "center_width") == 0) {
            double val;
			if (snd_config_get_real(n, &val) < 0) {
				SNDERR("Invalid type for %s", id);
				return -EINVAL;
			}
			center_width = (float)val;  // center_width [0..1] distributes the center information towards the front left/right channels, 1=full distribution, 0=no distribution
			if (center_width < 0.0 || center_width > 1.0) {
				SNDERR("center_width must be between 0.0 and 1.0");
				return -EINVAL;
			}
			continue;
		}
         if (strcmp(id, "dimension") == 0) {
            double val;
			if (snd_config_get_real(n, &val) < 0) {
				SNDERR("Invalid type for %s", id);
				return -EINVAL;
			}
			dimension = (float)val;  //  dimension [0..1] moves the soundfield backwards, 0=front, 1=side
			if (dimension < -0.5 || dimension > 1.0) {
				SNDERR("dimension must be between -0.5 and 1.0");
				return -EINVAL;
			}
			continue;
		}       
         if (strcmp(id, "adaption_rate") == 0) {
            double val;
			if (snd_config_get_real(n, &val) < 0) {
				SNDERR("Invalid type for %s", id);
				return -EINVAL;
			}
			adaption_rate = (float)val;  // adaption_rate [0..1] determines how fast the steering gets adapted, 1=instantaneous, 0.1 = very slow adaption
			if (adaption_rate < 0.0 || adaption_rate > 1.0) {
				SNDERR("adaption_rate must be between 0.0 and 1.0");
				return -EINVAL;
			}
			continue;
		}   
        if (strcmp(id, "phase_mode") == 0) {
            long val;
			if (snd_config_get_integer(n, &val) < 0) {
				SNDERR("Invalid type for %s", id);
				return -EINVAL;
			}
			phase_mode = (unsigned int)val;  
			if (phase_mode < 0 || phase_mode > 3) {
				SNDERR("phase_mode must be between 0 and 3");
				return -EINVAL;
			}
			continue;
		}   
        if (strcmp(id, "linear_steering") == 0) {
            long val;
			if (snd_config_get_integer(n, &val) < 0) {
				SNDERR("Invalid type for %s", id);
				return -EINVAL;
			}
			linear_steering = (unsigned int)val;  // 
			if (linear_steering != 0 && linear_steering != 1) {
				SNDERR("linear_steering must be 0 or 1");
				return -EINVAL;
			}
			continue;
		}   
        if (strcmp(id, "front_separation") == 0) {
            double val;
			if (snd_config_get_real(n, &val) < 0) {
				SNDERR("Invalid type for %s", id);
				return -EINVAL;
			}
			front_separation = (float)val;  // front_separation [0..1.5]
			if (front_separation < 0.0 || front_separation > 1.5) {
				SNDERR("front_separation must be between 0.0 and 1.5");
				return -EINVAL;
			}
			continue;
		}        
        if (strcmp(id, "rear_separation") == 0) {
            double val;
			if (snd_config_get_real(n, &val) < 0) {
				SNDERR("Invalid type for %s", id);
				return -EINVAL;
			}
			rear_separation = (float)val;  // rear_separation [0..1.5]
			if (rear_separation < 0.0 || rear_separation > 1.5) {
				SNDERR("rear_separation must be between 0.0 and 1.5");
				return -EINVAL;
			}
			continue;
		}        
		SNDERR("Unknown field %s", id);
		return -EINVAL;
	}
    
    if (! sconf) {
		SNDERR("No slave configuration for freesurround pcm");
		return -EINVAL;
	}

	fsmix = calloc(1, sizeof(*fsmix));
	if (! fsmix) {
		SNDERR("cannot allocate");
		return -ENOMEM;
	}

	fsmix->rate = rate;
	fsmix->debug = debug;
	
	fsmix->bitrate = bitrate;
	fsmix->channels = channels;
	fsmix->format = format;
    
    fsmix->center_width = center_width;
    fsmix->dimension = dimension;
    fsmix->adaption_rate = adaption_rate;
    fsmix->phase_mode = phase_mode;
    fsmix->linear_steering = linear_steering;
    fsmix->front_separation = front_separation;
    fsmix->rear_separation = rear_separation;


	fsmix->ext.version = SND_PCM_IOPLUG_VERSION;
	fsmix->ext.name = "FreeSurround upmix plugin";
	fsmix->ext.callback = &fs_callback;
	fsmix->ext.private_data = fsmix;
	if (fsmix->debug)
		printf("playing through %s\n",fsmix->ext.name);
	err = snd_pcm_extplug_create(&fsmix->ext, name, root, sconf, stream, mode);
	if (err < 0) {
        free(fsmix);
		return err;
    }

	snd_pcm_extplug_set_param_minmax(&fsmix->ext,
					 SND_PCM_EXTPLUG_HW_CHANNELS,
					 1, 6);
	if (channels)
		snd_pcm_extplug_set_slave_param_minmax(&fsmix->ext,
						       SND_PCM_EXTPLUG_HW_CHANNELS,
						       channels, channels);
	else
		snd_pcm_extplug_set_slave_param_list(&fsmix->ext,
						     SND_PCM_EXTPLUG_HW_CHANNELS,
						     2, chlist);
	snd_pcm_extplug_set_param(&fsmix->ext, SND_PCM_EXTPLUG_HW_FORMAT,
				  SND_PCM_FORMAT_S16);
	snd_pcm_extplug_set_slave_param(&fsmix->ext, SND_PCM_EXTPLUG_HW_FORMAT,
					SND_PCM_FORMAT_S16);

	*pcmp = fsmix->ext.pcm;
	return 0;

}

SND_PCM_PLUGIN_SYMBOL(freesurround);


