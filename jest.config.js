/* eslint-disable @typescript-eslint/no-var-requires */
module.exports = {
    preset: 'ts-jest',
    testEnvironment: 'node',
    verbose: true,
    globals: {
        'ts-jest': {
            compiler: 'typescript',
            tsconfig: 'tsconfig-jest.json',
            diagnostics: {
                ignoreCodes: [151001],
            },
        },
    },
};
