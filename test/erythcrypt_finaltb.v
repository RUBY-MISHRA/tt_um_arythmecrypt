`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date: 08/24/2024 09:23:29 AM
// Design Name: 
// Module Name: erythcrypt_finaltb
// Project Name: 
// Target Devices: 
// Tool Versions: 
// Description: 
// 
// Dependencies: 
// 
// Revision:
// Revision 0.01 - File Created
// Additional Comments:
// 
//////////////////////////////////////////////////////////////////////////////////


module erythcrypt_finaltb();
    reg [7:0] I1, I2;
    reg [3:0] Control;
    reg CLK, Reset;
    wire [7:0] OUTPUT;

    erythcrypt_final ef1 (
        .I1(I1),
        .I2(I2),
        .CLK(CLK),
        .Reset(Reset),
        .Control(Control),
        .OUTPUT(OUTPUT)
    );

    initial begin
        // Initialize signals
        I1 = 8'b00000000;
        I2 = 8'b00000000;
        Reset = 1;
        CLK = 0;
        Control = 4'b0000;
        
        // Release reset after some time
        #20 Reset = 0;

        // Apply test vectors
        #100 Control = 4'b0001;
        I1 = 8'b00011110; // 30
        I2 = 8'b01000110; // 70
//        I1 = 30; // 20
//                I2 = 70; // 30

        #1000 Control = 4'b0010;
        I1 = 8'b10000000; // 128
        I2 = 8'b01000001; // 65

        #1000 Control = 4'b0011;
        I1 = 8'b00001111; // 15
        I2 = 8'b00000001; // 1

        #1000 Control = 4'b0100;
        I1 = 8'b00001111; // 15
        I2 = 8'b00000001; // 1

        #1000 Control = 4'b0101;
        I1 = 8'b00010110; // 22
        I2 = 8'b00001010; // 10

        #1000 Control = 4'b0110;
        I1 = 8'b10000111; // 135
        I2 = 8'b01111011; // 123

        #1000 Control = 4'b0111;
        I1 = 8'b00000101; // 5
        I2 = 8'b00000010; // 2

        #1000 Control = 4'b1000;
        I1 = 8'b00000010; // 2
        I2 = 8'b00001100; // 12

        #1000 Control = 4'b1001;
        I1 = 8'b00000101; // 5
        I2 = 8'b00001100; // 12

        #400 Control = 4'b1010;
        I1 = 8'b00011010; // 26
        I2 = 8'b10000000; // 128
//I1=32;
//I2=80;
        #500 Control = 4'b1011;
        I1 = 8'b00000010; // 2
        I2 = 8'b00000101; // 5

        // Stop the simulation after running all tests
        #1000 $stop;
    end

    // Clock generation
    always #10 CLK = ~CLK;

endmodule
