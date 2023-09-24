from manim import *
import os

sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..')))

import properties.data_loader as dataloader
import csv

# Video configuration
config = dataloader.config
output_csv = os.path.join(os.path.dirname(__file__), config["meta"]["output_file"])
frame_rate = 30


### End of video configuration ###

def csv_to_dict(csv_path):
    with open(csv_path, 'r') as file:
        reader = csv.reader(file)
        headers = next(reader) # Read the header row
        
        data = []
        for row in reader:
            record = {}
            for i, value in enumerate(row):
                record[headers[i]] = value
            data.append(record)
        return data
    raise "Unable to open " + csv_path

class SimulationOutput(ThreeDScene):
    def construct(self):
        frame_time = 1/frame_rate
        cur_timer = 0
        total_time = 0
        file_data = csv_to_dict(output_csv)
        cur_row_num = 0

        max_frames = 40
        cur_frame = 0

        sim_time_text = Text(f"sim time:").scale(0.33).move_to([-5.75,-3,0])
        sim_time_num = DecimalNumber(0).scale(0.33).move_to(sim_time_text.get_right() + [0.25,0,0])
        self.add(sim_time_text)
        self.add(sim_time_num)

        self.add(NumberPlane().set_opacity(0.2))

        while(cur_row_num < len(file_data) and cur_frame <= max_frames):
            if(cur_row_num > 0):
                dt = float(file_data[cur_row_num]['time']) - float(file_data[cur_row_num-1]['time'])
                total_time += dt
                cur_timer += dt

            if(cur_timer >= frame_time):
                # We have passed the threshold, render this frame.
                sim_time_num.set_value(total_time)

                cur_timer = cur_timer - frame_time
                cur_frame += 1
                self.wait(frame_time)

            cur_row_num += 1
        
