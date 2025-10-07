//Include the Butterworth designer
#include <butterworth_design.hpp>

//Only bandpass filter
using FilterDesigner::ButterworthBandpass;

//Filter parameters
const float g_samples = 200.0;  //Sample rate frequency [Hz]
const float g_lowcut  = 15.0;   //Low cut frequency     [Hz]
const float g_highcut = 25.0;   //High cut frequency    [Hz]
const short g_order   = 4;      //Filter order

//DC offset for plot visualization
const float g_input_offset  = 3.f;
const float g_output_offset = 6.5f;

//Create a Butterworth band-pass design.
ButterworthBandpass<float, g_order> filter(g_samples, g_lowcut, g_highcut);

void setup() {
  //Serial plotter for debug
  Serial.begin(115200);

  //If filter synthesis fail, do not proceed
  if (!filter.designSuccess()) while (true) {}

}

void loop() {
  //Time control variables
  static float interval_micros = 0;
  static unsigned long last_micros = 0;
  
  //Sample rample control
  if (micros() - last_micros >= (1000000.0 / g_samples)) {
    last_micros = micros();

    //Sinusoidal signals
    float s1 = sin(1.0 * PI * 10.0 * interval_micros);  // 10 Hz
    float s2 = sin(1.0 * PI * 20.0 * interval_micros);  // 20 Hz
    float s3 = sin(1.0 * PI * 30.0 * interval_micros);  // 30 Hz

    //Merge signals
    float raw_signal = s1 + s2 + s3;

    //​Generates the output signal using a raw signal through biquad cascade.
    float out_signal = raw_signal;
    for (size_t index = 0; index < filter.m_num_sos; index++)
      out_signal = filter.m_sos_sections[index].process(out_signal);
    
    //Raw signal and output signal plot with some DC offset
    Serial.print("Input:");
    Serial.print(g_input_offset + raw_signal);
    Serial.print(" Output:");
    Serial.println(g_output_offset + out_signal);

    //Update time interval
    interval_micros += 1.f / g_samples;
  }
}