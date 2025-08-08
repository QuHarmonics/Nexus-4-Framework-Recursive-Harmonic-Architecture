#!/usr/bin/env python
#
# Copyright 2021 Sweet Tech AS
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

import random 
from PyInquirer import prompt
import time
import handy

def linearcmd_build(duration,position):
    # Add message to payload
    payload = handy.handyplug.Payload()
    message = payload.Messages.add()
    message.LinearCmd.Id = random.randint(0,4294967295)
    message.LinearCmd.DeviceIndex = 0
    vector = message.LinearCmd.Vectors.add()
    vector.Index = 0
    vector.Duration = duration
    vector.Position = position
    # Convert payload to string
    return payload.SerializeToString()

def inputvalidator(input,lowvalue,highvalue,outtype):
    try:
        return outtype(input)>=lowvalue and outtype(input)<=highvalue
    except ValueError:
        return "Invalid input"

def linearcmd_test(handyclient,duration,position):
    response = handyclient.send_data(linearcmd_build(duration,position))
    resp_payload = handy.handyplug.Payload()
    resp_payload.ParseFromString(response)
    resp_isok = False
    for resp_message in resp_payload.Messages:
        if resp_message.HasField("Ok"):
            resp_isok = True
            break
    if not resp_isok:
        raise RuntimeError(
            "Received incorrect response")

def main():

    # Create a new handy BLE client object
    print("Creating new Handy BLE client")
    handyclient = handy.BLEClient("The Handy")

    # Ask for commands
    while True:
        questions = [
            {
                'type': 'list',
                'name': 'command',
                'message': 'Choose a command',
                'choices': [
                    'Connect',
                    'Disconnect',
                    'LinearCmd',
                    'Sequential LinearCmd test'
                ]
            }
        ]
        answer = prompt(questions)

        if (answer['command'] == "Connect"):
            print("Attempting BLE connection")
            try:
                handyclient.connect()
            except RuntimeError as e:
                print(e)
        elif (answer['command'] == "Disconnect"):
            print("Disconnecting")
            try:
                handyclient.disconnect()
            except RuntimeError as e:
                print(e)
        elif (answer['command'] == "LinearCmd"):
            questions = [
                {
                    'type': 'input',
                    'name': 'duration',
                    'message': 'Duration? (ms)',
                    'validate': lambda val: inputvalidator(val,0,4294967295,int)
                },
                {
                    'type': 'input',
                    'name': 'position',
                    'message': 'Position? (0=0%, 1=100%)',
                    'validate': lambda val: inputvalidator(val,0,1,float)
                }
            ]
            answer = prompt(questions)
            try:
                linearcmd_test(handyclient,int(answer['duration']),float(answer['position']))
            except RuntimeError as e:
                print(e)
        elif (answer['command'] == "Sequential LinearCmd test"):
            try:
                linearcmd_test(handyclient,2000,1)
                time.sleep(1)
                linearcmd_test(handyclient,1000,0)
                time.sleep(0.5)
                linearcmd_test(handyclient,500,1)
                time.sleep(0.25)
                linearcmd_test(handyclient,250,0)
            except RuntimeError as e:
                print(e)

if __name__ == '__main__':
    main()# Re-import necessary libraries
import numpy as np

# Constants
planck_constant = 6.62607015e-34  # Joule·seconds
speed_of_light = 299792458  # meters per second
boltzmann_constant = 1.380649e-23  # Joules per Kelvin

# Cosmic Microwave Background (CMB) Peak Frequency
cmb_temperature = 2.725  # Kelvin
cmb_peak_frequency = (boltzmann_constant * cmb_temperature) / planck_constant

# Harmonic Resonance and Stability Calculation
harmonic_constant = 0.35  # Samson's Law factor

# Calculate the stabilization frequency based on recursive harmonic compression
cmb_harmonic_stabilization = cmb_peak_frequency * harmonic_constant

# Prepare results
cmb_peak_frequency, cmb_harmonic_stabilization
